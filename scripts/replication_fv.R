###############################################################
###############################################################
# Legacies of Resistance:
# Mobilization Against Organized Crime in Mexico
# 
# Javier Osorio
# 8/26/2020
#
#
# CONTENT:
# 1.  SET UP 
# 2.  PREAMBLE FOR SPATIAL ANALYSIS
# 3.  OLS MODEL
#     - Table 1 in the paper
#     - Table A3 in the Appendix
#     - Figure 3 in the paper
# 4.  2SLS MODEL
#     - Table 2 in the paper
#     - Table A5 in the Appendix
# 5.  SPATIAL 2SLS REGRESSION MODEL
#     - Table 3 in the paper 
#     - Table A6 in the Appendix 
# 6   APPENDIX
# 6.1 DESCRIPTIVE STATISTICS
#     - Table A1 in the Appendix
# 6.2 Logit model 
#     - Table A2 in the Appendix
# 6.3 SPATIAL REGRESSION MODEL
#     - Table A4 in the Appendix
# 6.4 Alternative explanations - inequality, crime, DTOs
#     - Table A7 in the Appendix
# 6.5 AIC 
#     - Table A8 in the Appendix
# 6.6 Alternative measure of Cristeros from Meyer 
#     - Table A9 in the Appendix
#     - Figure A2 in the Appendix
# 6.7 Robustness check for alternative explanations - using Brigades variable
# 6.8 Robustness check for alternative explanations - using Cristeros variable
#
###############################################################
###############################################################










###############################################################
###############################################################
# 1. SET UP 
###############################################################
###############################################################

# Load packages
library(foreign)  
library(sp) 
library(sphet) 
library(spdep)
library(AER)
library(maptools)
library(xtable)
library(stargazer)
library(texreg)
library(xtable)
library(dotwhisker)
library(broom)
library(dplyr)
library(sandwich)
library(lmtest)


##########################
# Set working directory
setwd("E:/Dropbox/Mexico_Vigilante_Project/Data/Autodefensas_reloaded/Data/Historical_Legacies_replication/data")

##########################
# Get the bishops data
bishops <- read.dta("./bishopsforR.dta") 










###############################################################
###############################################################
# 2. PREAMBLE FOR SPATIAL ANALYSIS
###############################################################
###############################################################


##########################
# Create binary weights spatial matrix

# Upload the data
muns <- readShapePoly("./muns/MUNICIPIOS.shp")

# Find contiguous neighborus and drop links
nb_q <- poly2nb(muns)
nb_q

# Check symetry
is.symmetric.nb(nb_q)

# Plot nodes
coords <- coordinates(muns)
plot(nb_q, coords, col="grey")


# Create spatial matrix as listw object - Binary weights
nb_B <- nb2listw(nb_q, style="B", zero.policy=TRUE)
nb_B$style
nb_B

# Create spatial matrix as listw object - Row-standardized weights
nb_R <- nb2listw(nb_q)
nb_R$style
nb_R


############################
#  Create inverse distance spatial matrix
dsts <-nbdists(nb_q,coordinates(coords))
idw<-lapply(dsts, function(x) 1/(x/1000))

# Create Inverse spatial matrix as listw object - Binary weights
nb_InvB <- nb2listw(nb_q, glist = idw, style="B", zero.policy=TRUE)
nb_InvB$style
nb_InvB

# Create Inverse spatial matrix as listw object - Row standardized weights
nb_InvR <- nb2listw(nb_q, glist = idw)
nb_InvR$style
nb_InvR










###############################################################
###############################################################
# 3. OLS MODEL
###############################################################
###############################################################

# Run different OLS model specifications by adding up covariates

# Model 1
ols1<-lm(selfdefense_presence2 ~ brig_all , 
         data=bishops) 
summary(ols1)
ols1r<-coeftest(ols1, vcov = vcovHC(ols1, "HC1"))    # robust SE as in Stata 
ols1r

# Model 2
ols2<-lm(selfdefense_presence2 ~ brig_all + 
           Railways + telegraphs + 
           rurales + morelos + mina + hidalgo + guerrero + french + 
           dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 , 
         data=bishops) 
summary(ols2)
ols2r<-coeftest(ols2, vcov = vcovHC(ols2, "HC1"))    # robust SE as in Stata 
ols2r

# Model 3
ols3<-lm(selfdefense_presence2 ~ brig_all + 
           Railways + telegraphs + 
           rurales + morelos + mina + hidalgo + guerrero + french + 
           dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
           pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929  , 
         data=bishops) 
summary(ols3)
ols3r<-coeftest(ols3, vcov = vcovHC(ols3, "HC1"))    # robust SE as in Stata 
ols3r

# Model 4
ols4<-lm(selfdefense_presence2 ~ brig_all + 
           Railways + telegraphs + 
           rurales + morelos + mina + hidalgo + guerrero + french + 
           dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
           pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
           loc_XVI + fran  + dom  + aug  + jes , 
         data=bishops) 
summary(ols4)
ols4r<-coeftest(ols4, vcov = vcovHC(ols4, "HC1"))    # robust SE as in Stata 
ols4r

# Model 5
ols5<-lm(selfdefense_presence2 ~ brig_all + 
           Railways + telegraphs + 
           rurales + morelos + mina + hidalgo + guerrero + french + 
           dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
           pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
           loc_XVI + fran  + dom  + aug  + jes +  
           arch_zone + triple_all + chichimeca , 
         data=bishops) 
summary(ols5)
ols5r<-coeftest(ols5, vcov = vcovHC(ols5, "HC1"))    # robust SE as in Stata 
ols5r


############################
# Generate tables of results


# This corresponds to Table 1 in the paper
# Abbreviated results
stargazer(ols1r, ols2r, ols3r, ols4r, ols5r, type="latex",   # models w robust SE
                    dep.var.labels=c("Autodefensas"),
          covariate.labels=c("Cristero Brigades"),
          keep = "brig_all",
          out = "../results/ols_brig.tex")


# This corresponds to Table A3 in the Appendix
# Full table of results
stargazer(ols1r, ols2r, ols3r, ols4r, ols5r, type="latex",        # models w robust SE    
          dep.var.labels=c("Autodefensas"),
          covariate.labels=c("Cristero Brigades",
                             "Railways", "Telegraphs",
                             "Rurales", "Morelos insurg.", "Mina insurg.", "Hidalgo insurg.", "Guerrero insurg.", "French invasion",  
                             "Distance to state capital", "Elevation", "Distance to railroad", "Gulf", "North", "Pacific",
                             "Pop. density 1921", "Rural pop. 1921", "Catholic pop. 1921", "Agricultural area 1928", "Federal employees 1929", "Police oficers 1929", "Iliteracy rate 1929", "Members per union 1929",
                             "Localities XVI Century", "Franciscan mission", "Dominican mission", "Augustine mission", "Jesuit mission", 
                             "Archeological zone","Tripple Alliance", "Chichimeca"),
          out = "../results/ols_all.tex" )



############################
# Generate coefficients plot 


# This corresponds to Figure 3 in the paper

# Extract Betas and SEs from each OLS model specification 
b1  <- ols1r[[2, 1,  exact = TRUE]]
se1 <- ols1r[[2, 2,  exact = TRUE]]
b2  <- ols2r[[2, 1,  exact = TRUE]]
se2 <- ols2r[[2, 2,  exact = TRUE]]
b3  <- ols3r[[2, 1,  exact = TRUE]]
se3 <- ols3r[[2, 2,  exact = TRUE]]
b4  <- ols4r[[2, 1,  exact = TRUE]]
se4 <- ols4r[[2, 2,  exact = TRUE]]
b5  <- ols5r[[2, 1,  exact = TRUE]]
se5 <- ols5r[[2, 2,  exact = TRUE]]

# Create a data frame of Betas and SEs
m<-cbind(c("Brigades", "Brigades","Brigades","Brigades","Brigades"),
         c(b1, b2, b3, b4, b5),
         c(se1, se2, se3, se4, se5),
         c("M1. No covariates", "M2: M1 + Infrastructure, \n \t\t Armed campaigns, \n \t\t Geography", "M3: M2 + Socio-Demographic", "M4: M3 + Colonial", "M5: M4 + Pre-Colonial"))

m<-data.frame(m)

names(m) <- c("Variable", "Coefficient", "SE", "Model")

# Specify the width of confidence intervals
interval1 <- -qnorm((1-0.9)/2)  # 90% multiplier
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# Make Coefficient and SE as numeric variables
m$Coefficient=as.numeric(m$Coefficient)
m$SE=as.numeric(m$SE)

# Plot
pdf("../graphs/ols_plot.pdf",  width=4, height=2) 
zp1 <- ggplot(m, aes(colour = Model))
zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
zp1 <- zp1 + geom_linerange(aes(x = Variable, ymin = Coefficient - SE*interval1,
                                ymax = Coefficient + SE*interval1),
                            lwd = 1, position = position_dodge(width = 1/2))
zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2),
                             lwd = 1/2, position = position_dodge(width = 1/2), # The trick to these is position_dodge().
                             shape = 21, fill = "WHITE")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + theme(axis.text.y=element_blank())
zp1 <- zp1 + labs(x = "Cristero brigades")
zp1
print(zp1) 
dev.off()










###############################################################
###############################################################
# 4. 2SLS MODEL
###############################################################
###############################################################

############################
# Instrumental Variables model

# First Stage
iv1_f <- lm(brig_all ~ bishop_dist100K +
              Railways + telegraphs + 
              rurales + morelos + mina + hidalgo + guerrero + french + 
              dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
              pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
              loc_XVI + fran  + dom  + aug  + jes +  
              arch_zone + triple_all + chichimeca , 
            data=bishops)
summary(iv1_f)
iv1_fr<-coeftest(iv1_f, vcov = vcovHC(iv1_f, "HC1"))    # robust SE as in Stata 
iv1_fr


# Second stage
iv1_s<-ivreg(selfdefense_presence2 ~ brig_all + 
               Railways + telegraphs + 
               rurales + morelos + mina + hidalgo + guerrero + french + 
               dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
               pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
               loc_XVI + fran  + dom  + aug  + jes +  
               arch_zone + triple_all + chichimeca 
             | Railways + telegraphs + 
               rurales + morelos + mina + hidalgo + guerrero + french + 
               dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
               pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
               loc_XVI + fran  + dom  + aug  + jes +  
               arch_zone + triple_all + chichimeca + bishop_dist100K , 
             data=bishops)  #, lag.instr=TRUE
summary(iv1_s, vcov = sandwich, diagnostics = TRUE)
iv1_sr<-coeftest(iv1_s, vcov = vcovHC(iv1_s, "HC1"))    # robust SE as in Stata 
iv1_sr

# Reduced Form
iv1_rf <- lm(selfdefense_presence2 ~ bishop_dist100K +
               Railways + telegraphs + 
               rurales + morelos + mina + hidalgo + guerrero + french + 
               dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
               pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
               loc_XVI + fran  + dom  + aug  + jes +  
               arch_zone + triple_all + chichimeca , 
             data=bishops)
summary(iv1_rf)
iv1_rfr<-coeftest(iv1_rf, vcov = vcovHC(iv1_rf, "HC1"))    # robust SE as in Stata 
iv1_rfr



############################
# Generate tables of results


# This corresponds to Table 2 in the paper
# Abbreviated results
stargazer(iv1_fr, iv1_sr, iv1_rfr , type="latex",                               # Results w robust SE
          dep.var.labels=c("Brigades","Autodefensas","Autodefensas"),
          covariate.labels=c("Distance to Pro-Cristero Bishops", "Cristero Brigades","Distance to Pro-Cristero Bishops"),
          keep = c("bishop_dist100K","brig_all" , "bishop_dist100K","Intercept"),
          out = "../results/iv_brig.tex")


# This corresponds to Table A5 in the Appendix
# Full table of results
stargazer(iv1_fr, iv1_sr, iv1_rfr , type="latex",                               # Results w robust SE
          dep.var.labels=c("Brigades","Autodefensas","Autodefensas"),
          covariate.labels=c("Distance to Pro-Cristero Bishops","Cristero Brigades",
                             "Railways", "Telegraphs",
                             "Rurales", "Morelos insurg.", "Mina insurg.", "Hidalgo insurg.", "Guerrero insurg.", "French invasion",  
                             "Distance to state capital", "Elevation", "Distance to railroad", "Gulf", "North", "Pacific",
                             "Pop. density 1921", "Rural pop. 1921", "Catholic pop. 1921", "Agricultural area 1928", "Federal employees 1929", "Police oficers 1929", "Iliteracy rate 1929", "Members per union 1929",
                             "Localities XVI Century", "Franciscan mission", "Dominican mission", "Augustine mission", "Jesuit mission", 
                             "Archeological zone","Tripple Alliance", "Chichimeca"),
          out = "../results/iv_brig_all.tex")










##################################################################
###############################################################
# 5. SPATIAL 2SLS REGRESSION MODEL
###############################################################
###############################################################




########################################################
# Run Spatial IV Model (siv)
# Full model - contiguous inverse binary spatial matrix
# Model = spatial Lag
# Heteroskedasticity


############################
# First stage
siv1_fh<-spreg(brig_all ~ bishop_dist100K +
                 Railways + telegraphs + 
                 rurales + morelos + mina + hidalgo + guerrero + french + 
                 dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                 pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                 loc_XVI + fran  + dom  + aug  + jes +  
                 arch_zone + triple_all + chichimeca , 
               data=bishops, listw=nb_InvB, endog = NULL, instruments = NULL, 
               lag.instr = TRUE, initial.value = 0.2, model = "lag", het=TRUE)
summary(siv1_fh)

sink("../results/siv1_fh.txt")
summary(siv1_fh)
sink() 




############################
# Second stage
siv1_sh<-spreg(selfdefense_presence2 ~  Railways + telegraphs + 
                 rurales + morelos + mina + hidalgo + guerrero + french + 
                 dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                 pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                 loc_XVI + fran  + dom  + aug  + jes +  
                 arch_zone + triple_all + chichimeca , 
               data=bishops, listw=nb_InvB, endog = ~ brig_all, instruments = ~ bishop_dist100K, 
               initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_sh)

# export
sink("../results/siv1_sh.txt")
summary(siv1_sh)
sink() 





############################
# Reduced form
siv1_rfh<-spreg(selfdefense_presence2 ~ bishop_dist100K + 
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca , 
                data=bishops, listw=nb_InvB,  
                initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_rfh)

# export
sink("../results/siv1_rfh.txt")
summary(siv1_rfh)
sink() 




############################
# Tables of results
#
# Unfortunately, the "stargazer" package does not recognize the "sphet" object class. 
# So, we cannot automatically generate a LaTeX table of results.
# Instead, we compiled the results of each model (siv1_fh, siv1_sh, siv1_rfh)
# into an Excel file, "siv_results.xlsx" to manually generate 
# the table of results for LaTeX.
# 
# This file is used to generate:
# - Table 3 in the paper (abbreviated results) and 
# - Table A6 in the Appendix (full results)










###############################################################
###############################################################
# 6 APPENDIX
###############################################################
###############################################################





###############################################################
# 6.1. DESCRIPTIVE STATISTICS
###############################################################


# Replicate Table A1 in the Appendix 

# Generate data for descriptive statistics
desc <- subset(bishops, select= c("selfdefense_presence2" ,   
                                  "brig_all", "cristeros",
                                  "bishop_dist100K",
                                  "Railways" , "telegraphs" , 
                                  "dist_state_capital_log" ,"elevation_log" , "dist_rail_k" , "gulf_3" , "North_3" , "Pacific_3" , 
                                  "rurales" , "morelos" , "mina" , "hidalgo" , "guerrero" , "french" , 
                                  "pop_den_1921" , "rural_p_1921" , "cathol_p_1921"  , "agri_area_1928" , "gov_fed_10K_1292" ,  "police_10K_1928" ,  "analfa_1921" ,  "mem_per_union_1929" , 
                                  "loc_XVI" , "fran"  , "dom"  , "aug"  , "jes" ,  
                                  "arch_zone" , "triple_all" , "chichimeca"))
desc <-data.frame(desc)

# Print summary statistics
summary(desc)

# Generate Table 4 reporting the summary statistics
stargazer(desc, type = "latex", title="Descriptive statistics", 
          summary.stat = c("mean", "sd", "min", "max", "n"),
          digits=1, out="../results/desc.tex",
          covariate.labels=c("Autodefensas", 
                             "Cristero Brigades", "Cristeros",
                             "Distance to pro-Cristero Bishops",
                             "Railways", "Telegraphs",
                             "Distance to state capital", "Elevation", "Distance to railroad", "Gulf", "North", "Pacific",
                             "Rurales", "Morelos insurg.", "Mina insurg.", "Hidalgo insurg.", "Guerrero insurg.", "French invasion",  
                             "Pop. density 1921", "Rural pop. 1921", "Catholic pop. 1921", "Agricultural area 1928", "Federal employees 1929", "Police oficers 1929", "Iliteracy rate 1929", "Members per union 1929",
                             "Localities XVI Century", "Franciscan mission", "Dominican mission", "Augustine mission", "Jesuit mission", 
                             "Archeological zone", "Tripple Alliance", "Chichimeca")
)





###############################################################
# 6.2 Logit model 
###############################################################


# Replicates Table A2 in the Appendix

# Function for calculating robust SE as in Stata
sandwich1 <- function(object, ...) sandwich(object) * nobs(object) / (nobs(object) - 1)


############################
# Run logit models

# Logit Model 1
logit1<-glm(selfdefense_presence2 ~ brig_all , 
            data=bishops,
            family=binomial(link="logit")) 
summary(logit1)
logit1r <- coeftest(logit1, vcov = sandwich1)    # Calculate robust SE
logit1r

# Logit Model 2
logit2<-glm(selfdefense_presence2 ~ brig_all + 
              Railways + telegraphs + 
              rurales + morelos + mina + hidalgo + guerrero + french + 
              dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 , 
            data=bishops,
            family=binomial(link="logit")) 
summary(logit2)
logit2r <- coeftest(logit2, vcov = sandwich1)    # Calculate robust SE
logit2r

# Logit Model 3
logit3<-glm(selfdefense_presence2 ~ brig_all + 
              Railways + telegraphs + 
              rurales + morelos + mina + hidalgo + guerrero + french + 
              dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
              pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929  , 
            data=bishops,
            family=binomial(link="logit")) 
summary(logit3)
logit3r <- coeftest(logit3, vcov = sandwich1)    # Calculate robust SE
logit3r

# Logit Model 4
logit4<-glm(selfdefense_presence2 ~ brig_all + 
              Railways + telegraphs + 
              rurales + morelos + mina + hidalgo + guerrero + french + 
              dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
              pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
              loc_XVI + fran  + dom  + aug  + jes , 
            data=bishops,
            family=binomial(link="logit")) 
summary(logit4)
logit4r <- coeftest(logit4, vcov = sandwich1)    # Calculate robust SE
logit4r

# Logit Model 5
logit5<- glm(selfdefense_presence2 ~ brig_all + 
               Railways + telegraphs + 
               rurales + morelos + mina + hidalgo + guerrero + french + 
               dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
               pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
               loc_XVI + fran  + dom  + aug  + jes +  
               arch_zone + triple_all + chichimeca , 
             data=bishops,
             family=binomial(link="logit")) 
summary(logit5)
logit5r <- coeftest(logit5, vcov = sandwich1)    # Calculate robust SE
logit5r




############################
# Generate table of results


# Replicate Table A2 in the Appendix

stargazer(logit1r, logit2r, logit3r, logit4r, logit5r, type="latex",    # w robust SE          dep.var.labels=c("Autodefensas"),
          covariate.labels=c("Cristero Brigades",
                             "Railways", "Telegraphs",
                             "Rurales", "Morelos insurg.", "Mina insurg.", "Hidalgo insurg.", "Guerrero insurg.", "French invasion",  
                             "Distance to state capital", "Elevation", "Distance to railroad", "Gulf", "North", "Pacific",
                             "Pop. density 1921", "Rural pop. 1921", "Catholic pop. 1921", "Agricultural area 1928", "Federal employees 1929", "Police oficers 1929", "Iliteracy rate 1929", "Members per union 1929",
                             "Localities XVI Century", "Franciscan mission", "Dominican mission", "Augustine mission", "Jesuit mission", 
                             "Archeological zone","Tripple Alliance", "Chichimeca"),
          out = "../results/logit_all.tex" )










###############################################################
# 6.3 SPATIAL REGRESSION MODEL
###############################################################


# Replicate Table A4 in the Appendix


############################
# Spatial regression
# The "het" option already considers hetoroskedasticity


# Model 1
s1<-spreg(selfdefense_presence2 ~ brig_all  ,
          data=bishops, listw=nb_InvB,  
          initial.value = 0.2, model = "lag", het=TRUE)  
summary(s1)

sink("../results/s1.txt")
summary(s1)
sink() 


# Model 2
s2<-spreg(selfdefense_presence2 ~ brig_all + 
            Railways + telegraphs + 
            rurales + morelos + mina + hidalgo + guerrero + french + 
            dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 ,
          data=bishops, listw=nb_InvB,  
          initial.value = 0.2, model = "lag", het=TRUE)  
summary(s2)

sink("../results/s2.txt")
summary(s2)
sink() 


# Model 3
s3<-spreg(selfdefense_presence2 ~ brig_all + 
            Railways + telegraphs + 
            rurales + morelos + mina + hidalgo + guerrero + french + 
            dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
            pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 ,
          data=bishops, listw=nb_InvB,  
          initial.value = 0.2, model = "lag", het=TRUE)  
summary(s3)

sink("../results/s3.txt")
summary(s3)
sink() 


# Model 4
s4<-spreg(selfdefense_presence2 ~ brig_all + 
            Railways + telegraphs + 
            rurales + morelos + mina + hidalgo + guerrero + french + 
            dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
            pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
            loc_XVI + fran  + dom  + aug  + jes,
          data=bishops, listw=nb_InvB,  
          initial.value = 0.2, model = "lag", het=TRUE)  
summary(s4)

sink("../results/s4.txt")
summary(s4)
sink() 


# Model 5
s5<-spreg(selfdefense_presence2 ~ brig_all + 
            Railways + telegraphs + 
            rurales + morelos + mina + hidalgo + guerrero + french + 
            dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
            pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
            loc_XVI + fran  + dom  + aug  + jes +  
            arch_zone + triple_all + chichimeca,
          data=bishops, listw=nb_InvB,  
          initial.value = 0.2, model = "lag", het=TRUE)  
summary(s5)

sink("../results/s5.txt")
summary(s5)
sink() 


############################
# Tables of results
#
# Unfortunately, the "stargazer" package does not recognize the "sphet" object class. 
# So, we cannot automatically generate a LaTeX table of results.
# Instead, we compiled the results of each model (s1, s2, s3, s4, s5)
# into an Excel file, "s_results.xlsx" to manually generate 
# the table of results for LaTeX.
# 
# This file is used to generate Table A4 in the Appendix (full results)










###############################################################
# 6.4 Alternative explanations - inequality, crime, DTOs 
###############################################################

# Replicate Table A7 in the Appendix


##########################
# Get the bishops data
inequality <- read.dta("./phillipsvigilantesCPS.dta") 


# Merge databases
bishops2 <- merge(bishops, inequality, by.x="cve_mun_old", by.y="idedomun")


##########################
# Model Inequality (full specification)
ineq1<-lm( gini  ~ brig_all + 
             Railways + telegraphs + 
             rurales + morelos + mina + hidalgo + guerrero.x + french + 
             dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
             pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
             loc_XVI + fran  + dom  + aug  + jes +  
             arch_zone + triple_all + chichimeca , 
           data=bishops2) 
summary(ineq1)
ineq1r<-coeftest(ineq1, vcov = vcovHC(ineq1, "HC1"))    # robust SE as in Stata 
ineq1r


##########################
# Model Crime -  (full specification)
crime1<-lm( totalao2013  ~ brig_all + 
              Railways + telegraphs + 
              rurales + morelos + mina + hidalgo + guerrero.x + french + 
              dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
              pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
              loc_XVI + fran  + dom  + aug  + jes +  
              arch_zone + triple_all + chichimeca , 
            data=bishops2) 
summary(crime1)
crime1r<-coeftest(crime1, vcov = vcovHC(crime1, "HC1"))    # robust SE as in Stata 
crime1r


##########################
# Model DTOs (full specification) 
dtos1<-lm( total_all_dto10  ~ brig_all + 
             Railways + telegraphs + 
             rurales + morelos + mina + hidalgo + guerrero.x + french + 
             dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
             pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
             loc_XVI + fran  + dom  + aug  + jes +  
             arch_zone + triple_all + chichimeca , 
           data=bishops2) 
summary(dtos1)
dtos1r<-coeftest(dtos1, vcov = vcovHC(dtos1, "HC1"))    # robust SE as in Stata 
dtos1r


############################
# Generate  Table A7 in the Appendix

stargazer(ineq1r, crime1r, dtos1r, type="latex",        # models w robust SE    
          covariate.labels=c("Cristero Brigades",
                             "Railways", "Telegraphs",
                             "Rurales", "Morelos insurg.", "Mina insurg.", "Hidalgo insurg.", "Guerrero insurg.", "French invasion",  
                             "Distance to state capital", "Elevation", "Distance to railroad", "Gulf", "North", "Pacific",
                             "Pop. density 1921", "Rural pop. 1921", "Catholic pop. 1921", "Agricultural area 1928", "Federal employees 1929", "Police oficers 1929", "Iliteracy rate 1929", "Members per union 1929",
                             "Localities XVI Century", "Franciscan mission", "Dominican mission", "Augustine mission", "Jesuit mission", 
                             "Archeological zone","Tripple Alliance", "Chichimeca"),
          out = "../results/alternative.tex" )










###############################################################
# 6.5 AIC 
###############################################################

# Replicate Table A8 in the Appendix

############################################
# Calculate AIC statistic for the OLS 5 model
Sum1 <- summary(ols5)
RSS1 <- sum(Sum1$residuals^2)
K1 <- length(coef(ols5))-1
N1 <- length(ols5$residuals)
n1 <- N1 - K1 - ols5$df.residual

AIC1 = log(RSS1/n1) + (2*K1)/n1
AIC1 <-format(round(AIC1, 2), nsmall = 2)



############################################
# Calculate AIC statistic for the Spatial model 
Sum2 <- summary(s5)
RSS2 <- sum(Sum2$residuals^2)
K2 <- length(coef(s5))-1
N2 <- length(s5$residuals)
tdf2 <- N2-1                       # Total DF
mdf2 <- length(coef(s5))-1    # Model DF
rdf2 <- tdf2-mdf2                  # Residual DF
n2 <- N2 - K2 - rdf2

AIC2 = log(RSS2/n2) + (2*K2)/n2
AIC2 <-format(round(AIC2, 2), nsmall = 2)



############################################
# Calculate AIC statistic for the IV model - second stage
Sum3 <- summary(iv1_s)
RSS3 <- sum(Sum3$residuals^2)
K3 <- length(coef(iv1_s))-1
N3 <- length(iv1_s$residuals)
n3 <- N3 - K3 - iv1_s$df.residual

AIC3 = log(RSS3/n3) + (2*K3)/n3
AIC3 <-format(round(AIC3, 2), nsmall = 2)



############################################
# Calculate AIC statistic for the Spatial IV model - second stage
Sum4 <- summary(siv1_sh)
RSS4 <- sum(Sum4$residuals^2)
K4 <- length(coef(siv1_sh))-1
N4 <- length(siv1_sh$residuals)
tdf4 <- N4-1                       # Total DF
mdf4 <- length(coef(siv1_sh))-1    # Model DF
rdf4 <- tdf4-mdf4                  # Residual DF
n4 <- N4 - K4 - rdf4

AIC4 = log(RSS4/n4) + (2*K4)/n4
AIC4 <-format(round(AIC4, 2), nsmall = 2)




############################################
# Report AIC table

# Replicate Table A8 in the Appendix

header <-c("OLS", "Spatial", "IV 2nd stage", "Spatial IV 2nd stage")
AIC <- c(AIC1,AIC2,AIC3,AIC4)
AIC.table <-rbind(header, AIC)
AIC.table

# Generate LaTeX table
library(xtable)
xtable(AIC.table)










###############################################################
# 6.6 Alternative measure of Cristeros from Meyer 
###############################################################

# Replicate Table A9 in the Appendix

# Model 1
olsb1<-lm(selfdefense_presence2 ~ cristeros , 
          data=bishops) 
summary(olsb1)
olsb1r<-coeftest(olsb1, vcov = vcovHC(olsb1, "HC1"))    # robust SE as in Stata 
olsb1r

# Model 2
olsb2<-lm(selfdefense_presence2 ~ cristeros + 
            Railways + telegraphs + 
            rurales + morelos + mina + hidalgo + guerrero + french + 
            dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 , 
          data=bishops) 
summary(olsb2)
olsb2r<-coeftest(olsb2, vcov = vcovHC(olsb2, "HC1"))    # robust SE as in Stata 
olsb2r

# Model 3
olsb3<-lm(selfdefense_presence2 ~ cristeros + 
            Railways + telegraphs + 
            rurales + morelos + mina + hidalgo + guerrero + french + 
            dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
            pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929  , 
          data=bishops) 
summary(olsb3)
olsb3r<-coeftest(olsb3, vcov = vcovHC(olsb3, "HC1"))    # robust SE as in Stata 
olsb3r

# Model 4
olsb4<-lm(selfdefense_presence2 ~ cristeros + 
            Railways + telegraphs + 
            rurales + morelos + mina + hidalgo + guerrero + french + 
            dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
            pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
            loc_XVI + fran  + dom  + aug  + jes , 
          data=bishops) 
summary(olsb4)
olsb4r<-coeftest(olsb4, vcov = vcovHC(olsb4, "HC1"))    # robust SE as in Stata 
olsb4r

# Model 5
olsb5<-lm(selfdefense_presence2 ~ cristeros + 
            Railways + telegraphs + 
            rurales + morelos + mina + hidalgo + guerrero + french + 
            dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
            pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
            loc_XVI + fran  + dom  + aug  + jes +  
            arch_zone + triple_all + chichimeca , 
          data=bishops) 

summary(olsb5)
olsb5r<-coeftest(olsb5, vcov = vcovHC(olsb5, "HC1"))    # robust SE as in Stata 
olsb5r




############################
# Generate table of results

# Replicate Table A9 in the Appendix

# All coefficients
stargazer(olsb1r, olsb2r, olsb3r, olsb4r, olsb5r, type="latex",    # w robust SE
          dep.var.labels=c("Autodefensas"),
          covariate.labels=c("Cristeros",
                             "Railways", "Telegraphs",
                             "Rurales", "Morelos insurg.", "Mina insurg.", "Hidalgo insurg.", "Guerrero insurg.", "French invasion",  
                             "Distance to state capital", "Elevation", "Distance to railroad", "Gulf", "North", "Pacific",
                             "Pop. density 1921", "Rural pop. 1921", "Catholic pop. 1921", "Agricultural area 1928", "Federal employees 1929", "Police oficers 1929", "Iliteracy rate 1929", "Members per union 1929",
                             "Localities XVI Century", "Franciscan mission", "Dominican mission", "Augustine mission", "Jesuit mission", 
                             "Archeological zone","Tripple Alliance", "Chichimeca"),
          out = "../results/ols_cris_all.tex" )




############################
# Generate coefficients plot

# Replicate Figure A2 in the Appendix

#olsb1r, olsb2r, olsb3r, olsb4r, olsb5r

# Extract betas and SEs from OLD with robuts SE
b1  <- olsb1r[[2, 1,  exact = TRUE]]
se1 <- olsb1r[[2, 2,  exact = TRUE]]
b2  <- olsb2r[[2, 1,  exact = TRUE]]
se2 <- olsb2r[[2, 2,  exact = TRUE]]
b3  <- olsb3r[[2, 1,  exact = TRUE]]
se3 <- olsb3r[[2, 2,  exact = TRUE]]
b4  <- olsb4r[[2, 1,  exact = TRUE]]
se4 <- olsb4r[[2, 2,  exact = TRUE]]
b5  <- olsb5r[[2, 1,  exact = TRUE]]
se5 <- olsb5r[[2, 2,  exact = TRUE]]

# Create a dataframe
m<-cbind(c("Brigades", "Brigades","Brigades","Brigades","Brigades"),
         c(b1, b2, b3, b4, b5),
         c(se1, se2, se3, se4, se5),
         c("M1. No covariates", "M2: M1 + Infrastructure, \n \t\t Armed campaigns, \n \t\t Geography", "M3: M2 + Socio-Demographic", "M4: M3 + Colonial", "M5: M4 + Pre-Colonial"))
m<-data.frame(m)
names(m) <- c("Variable", "Coefficient", "SE", "Model")

# Specify the width of your confidence intervals
interval1 <- -qnorm((1-0.9)/2)  # 90% multiplier
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# Make Coefficient and SE as numeric variables
m$Coefficient=as.numeric(m$Coefficient)
m$SE=as.numeric(m$SE)

# Plot
pdf("../graphs/ols_cris_plot.pdf",  width=4, height=3) 
zp2 <- ggplot(m, aes(colour = Model))
zp2 <- zp2 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
zp2 <- zp2 + geom_linerange(aes(x = Variable, ymin = Coefficient - SE*interval1,
                                ymax = Coefficient + SE*interval1),
                            lwd = 1, position = position_dodge(width = 1/2))
zp2 <- zp2 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2),
                             lwd = 1/2, position = position_dodge(width = 1/2), # The trick to these is position_dodge().
                             shape = 21, fill = "WHITE")
zp2 <- zp2 + coord_flip() + theme_bw()
zp2 <- zp2 + theme(axis.text.y=element_blank())
zp2 <- zp2 + labs(x = "Cristeros")
print(zp2) 
dev.off()










###############################################################
# 6.7 Robustness check for alternative explanations - using Brigades variable
# DV: inequality, crime, DTOs 
# These models use the Brigades variable
###############################################################


############################
# This section conducts a robustness check for alternative explanations 
# including inequality, crime, and DTOs. 
# These results are not included in the Appendix.
# We just provide a summary table on Table A10.



########################################################
# 6.7.1 2SLS MODEL

############################
# Inequality


# First Stage
iv1_fi_B <- lm(brig_all ~ bishop_dist100K +
                 Railways + telegraphs + 
                 rurales + morelos + mina + hidalgo + guerrero.x + french + 
                 dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                 pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                 loc_XVI + fran  + dom  + aug  + jes +  
                 arch_zone + triple_all + chichimeca , 
               data=bishops2)
summary(iv1_fi_B)
iv1_fi_Br<-coeftest(iv1_fi_B, vcov = vcovHC(iv1_fi_B, "HC1"))    # robust SE as in Stata 
iv1_fi_Br


# Second stage
iv1_si_A<-ivreg(gini ~ brig_all +  
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca 
                | Railways + telegraphs +  
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca + bishop_dist100K , 
                data=bishops2)  #, lag.instr=TRUE
summary(iv1_si_A, vcov = sandwich, diagnostics = TRUE)
iv1_si_Ar<-coeftest(iv1_si_A, vcov = vcovHC(iv1_si_A, "HC1"))    # robust SE as in Stata 
iv1_si_Ar


# Reduced Form
iv1_rfi_A <- lm(gini ~ bishop_dist100K +
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca , 
                data=bishops2)
summary(iv1_rfi_A)
iv1_rfi_Ar<-coeftest(iv1_rfi_A, vcov = vcovHC(iv1_rfi_A, "HC1"))    # robust SE as in Stata 
iv1_rfi_Ar



############################
# Crime

# First Stage
iv1_fc_A <- lm(brig_all ~ bishop_dist100K +
                 Railways + telegraphs + 
                 rurales + morelos + mina + hidalgo + guerrero.x + french + 
                 dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                 pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                 loc_XVI + fran  + dom  + aug  + jes +  
                 arch_zone + triple_all + chichimeca , 
               data=bishops2)
summary(iv1_fc_A)
iv1_fc_Ar<-coeftest(iv1_fc_A, vcov = vcovHC(iv1_fc_A, "HC1"))    # robust SE as in Stata 
iv1_fc_Ar


# Second stage
iv1_sc_A<-ivreg(totalao2013 ~ brig_all + 
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca 
                | Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca + bishop_dist100K , 
                data=bishops2)  #, lag.instr=TRUE
summary(iv1_sc_A, vcov = sandwich, diagnostics = TRUE)
iv1_sc_Ar<-coeftest(iv1_sc_A, vcov = vcovHC(iv1_sc_A, "HC1"))    # robust SE as in Stata 
iv1_sc_Ar


# Reduced Form
iv1_rfc_A <- lm(totalao2013 ~ bishop_dist100K +
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca , 
                data=bishops2)
summary(iv1_rfc_A)
iv1_rfc_Ar<-coeftest(iv1_rfc_A, vcov = vcovHC(iv1_rfc_A, "HC1"))    # robust SE as in Stata 
iv1_rfc_Ar



############################
# DTOs

# First Stage
iv1_fd_A <- lm(brig_all ~ bishop_dist100K +
                 Railways + telegraphs + 
                 rurales + morelos + mina + hidalgo + guerrero.x + french + 
                 dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                 pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                 loc_XVI + fran  + dom  + aug  + jes +  
                 arch_zone + triple_all + chichimeca , 
               data=bishops2)
summary(iv1_fd_A)
iv1_fd_Ar<-coeftest(iv1_fd_A, vcov = vcovHC(iv1_fd_A, "HC1"))    # robust SE as in Stata 
iv1_fd_Ar


# Second stage
iv1_sd_A<-ivreg(total_all_dto10 ~ brig_all + 
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca 
                | Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca + bishop_dist100K , 
                data=bishops2)  #, lag.instr=TRUE
summary(iv1_sd_A, vcov = sandwich, diagnostics = TRUE)
iv1_sd_Ar<-coeftest(iv1_sd_A, vcov = vcovHC(iv1_sd_A, "HC1"))    # robust SE as in Stata 
iv1_sd_Ar


# Reduced Form
iv1_rfd_A <- lm(total_all_dto10 ~ bishop_dist100K +
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca , 
                data=bishops2)
summary(iv1_rfd_A)
iv1_rfd_Ar<-coeftest(iv1_rfd_A, vcov = vcovHC(iv1_rfd_A, "HC1"))    # robust SE as in Stata 
iv1_rfd_Ar




########################################################
# 6.7.2 Spatial 2SLS MODEL



si_A<-spreg(gini2 ~ brig_all + 
              Railways + telegraphs + 
              rurales + morelos + mina + hidalgo + guerrero.x + french + 
              dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
              pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
              loc_XVI + fran  + dom  + aug  + jes +  
              arch_zone + triple_all + chichimeca,
            data=bishops2, listw=nb_InvB,  
            initial.value = 0.2, model = "lag", het=TRUE)  
summary(si_A)

sink("../results/si_A.txt")
summary(si_A)
sink() 


sc_A<-spreg(totalao2013 ~ brig_all + 
              Railways + telegraphs + 
              rurales + morelos + mina + hidalgo + guerrero.x + french + 
              dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
              pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
              loc_XVI + fran  + dom  + aug  + jes +  
              arch_zone + triple_all + chichimeca,
            data=bishops2, listw=nb_InvB,  
            initial.value = 0.2, model = "lag", het=TRUE)  
summary(sc_A)

sink("../results/sc_A.txt")
summary(sc_A)
sink() 


sd_A<-spreg(total_all_dto10 ~ brig_all + 
              Railways + telegraphs + 
              rurales + morelos + mina + hidalgo + guerrero.x + french + 
              dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
              pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
              loc_XVI + fran  + dom  + aug  + jes +  
              arch_zone + triple_all + chichimeca,
            data=bishops2, listw=nb_InvB,  
            initial.value = 0.2, model = "lag", het=TRUE)  
summary(sd_A)

sink("../results/sd_A.txt")
summary(sd_A)
sink() 





########################################################
# 6.7.3 Spatial 2SLS MODEL



############################
# Inequality


# First stage
siv1_fhi_A<-spreg(brig_all ~ bishop_dist100K +
                    Railways + telegraphs + 
                    rurales + morelos + mina + hidalgo + guerrero.x + french + 
                    dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                    pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                    loc_XVI + fran  + dom  + aug  + jes +  
                    arch_zone + triple_all + chichimeca , 
                  data=bishops2, listw=nb_InvB, endog = NULL, instruments = NULL, 
                  lag.instr = TRUE, initial.value = 0.2, model = "lag", het=TRUE)
summary(siv1_fhi_A)

sink("../results/siv1_fhi_A.txt")
summary(siv1_fhi_A)
sink() 



## Impute data for gini in two observations
ginimean <- mean(bishops2$gini, na.rm = TRUE)
bishops2$gini2<-bishops2$gini
bishops2$gini2[is.na(bishops2$gini2)] <- ginimean


# Second stage
siv1_shi_A<-spreg(gini2 ~  Railways + telegraphs + 
                    rurales + morelos + mina + hidalgo + guerrero.x + french + 
                    dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                    pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                    loc_XVI + fran  + dom  + aug  + jes +  
                    arch_zone + triple_all + chichimeca , 
                  data=bishops2, listw=nb_InvB, endog = ~ brig_all, instruments = ~ bishop_dist100K, 
                  initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_shi_A)

# export
sink("../results/siv1_shi_A.txt")
summary(siv1_shi_A)
sink() 


# Reduced form
siv1_rfhi_A<-spreg(gini2 ~ bishop_dist100K + 
                     Railways + telegraphs + 
                     rurales + morelos + mina + hidalgo + guerrero.x + french + 
                     dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                     pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                     loc_XVI + fran  + dom  + aug  + jes +  
                     arch_zone + triple_all + chichimeca , 
                   data=bishops2, listw=nb_InvB,  
                   initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_rfhi_A)

# export
sink("../results/siv1_rfhi_A.txt")
summary(siv1_rfhi_A)
sink() 



############################
# Crime

# First stage
siv1_fhc_A<-spreg(brig_all ~ bishop_dist100K +
                    Railways + telegraphs + 
                    rurales + morelos + mina + hidalgo + guerrero.x + french + 
                    dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                    pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                    loc_XVI + fran  + dom  + aug  + jes +  
                    arch_zone + triple_all + chichimeca , 
                  data=bishops2, listw=nb_InvB, endog = NULL, instruments = NULL, 
                  lag.instr = TRUE, initial.value = 0.2, model = "lag", het=TRUE)
summary(siv1_fhc_A)

sink("../results/siv1_fhc_A.txt")
summary(siv1_fhc_A)
sink() 


# Second stage
siv1_shc_A<-spreg(totalao2013 ~  Railways + telegraphs + 
                    rurales + morelos + mina + hidalgo + guerrero.x + french + 
                    dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                    pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                    loc_XVI + fran  + dom  + aug  + jes +  
                    arch_zone + triple_all + chichimeca , 
                  data=bishops2, listw=nb_InvB, endog = ~ brig_all, instruments = ~ bishop_dist100K, 
                  initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_shc_A)

# export
sink("../results/siv1_shc_A.txt")
summary(siv1_shc_A)
sink() 


# Reduced form
siv1_rfhc_A<-spreg(totalao2013 ~ bishop_dist100K + 
                     Railways + telegraphs + 
                     rurales + morelos + mina + hidalgo + guerrero.x + french + 
                     dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                     pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                     loc_XVI + fran  + dom  + aug  + jes +  
                     arch_zone + triple_all + chichimeca , 
                   data=bishops2, listw=nb_InvB,  
                   initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_rfhc_A)

# export
sink("../results/siv1_rfhc_A.txt")
summary(siv1_rfhc_A)
sink() 



############################
# DTOs

# First stage
siv1_fhd_A<-spreg(brig_all ~ bishop_dist100K +
                    Railways + telegraphs + 
                    rurales + morelos + mina + hidalgo + guerrero.x + french + 
                    dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                    pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                    loc_XVI + fran  + dom  + aug  + jes +  
                    arch_zone + triple_all + chichimeca , 
                  data=bishops2, listw=nb_InvB, endog = NULL, instruments = NULL, 
                  lag.instr = TRUE, initial.value = 0.2, model = "lag", het=TRUE)
summary(siv1_fhd_A)

sink("../results/siv1_fhd_A.txt")
summary(siv1_fhd_A)
sink() 


# Second stage
siv1_shd_A<-spreg(total_all_dto10 ~  Railways + telegraphs + 
                    rurales + morelos + mina + hidalgo + guerrero.x + french + 
                    dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                    pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                    loc_XVI + fran  + dom  + aug  + jes +  
                    arch_zone + triple_all + chichimeca , 
                  data=bishops2, listw=nb_InvB, endog = ~ brig_all, instruments = ~ bishop_dist100K, 
                  initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_shd_A)

# export
sink("../results/siv1_shd_A.txt")
summary(siv1_shd_A)
sink() 


# Reduced form
siv1_rfhd_A<-spreg(total_all_dto10 ~ bishop_dist100K + 
                     Railways + telegraphs + 
                     rurales + morelos + mina + hidalgo + guerrero.x + french + 
                     dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                     pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                     loc_XVI + fran  + dom  + aug  + jes +  
                     arch_zone + triple_all + chichimeca , 
                   data=bishops2, listw=nb_InvB,  
                   initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_rfhd_A)

# export
sink("../results/siv1_rfhd_A.txt")
summary(siv1_rfhd_A)
sink() 










###############################################################
# 6.8 Robustness check for alternative explanations - using Cristeros variable
# DV: inequality, crime, DTOs 
# Using alternative Cristeros variable
###############################################################


############################
# This section conducts a robustness check for alternative explanations 
# including inequality, crime, and DTOs. 
# These results are not included in the Appendix



########################################################
# 6.8.1 2SLS MODEL

############################
# Inequality


# First Stage
iv1_fi_C <- lm(cristeros ~ bishop_dist100K +
                 Railways + telegraphs + 
                 rurales + morelos + mina + hidalgo + guerrero.x + french + 
                 dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                 pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                 loc_XVI + fran  + dom  + aug  + jes +  
                 arch_zone + triple_all + chichimeca , 
               data=bishops2)
summary(iv1_fi_C)
iv1_fi_Cr<-coeftest(iv1_fi_C, vcov = vcovHC(iv1_fi_C, "HC1"))    # robust SE as in Stata 
iv1_fi_Cr


# Second stage
iv1_si_C<-ivreg(gini ~ cristeros +  
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca 
                | Railways + telegraphs +  
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca + bishop_dist100K , 
                data=bishops2)  #, lag.instr=TRUE
summary(iv1_si_C, vcov = sandwich, diagnostics = TRUE)
iv1_si_Cr<-coeftest(iv1_si_C, vcov = vcovHC(iv1_si_C, "HC1"))    # robust SE as in Stata 
iv1_si_Cr


# Reduced Form
iv1_rfi_C <- lm(gini ~ bishop_dist100K +
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca , 
                data=bishops2)
summary(iv1_rfi_C)
iv1_rfi_Cr<-coeftest(iv1_rfi_C, vcov = vcovHC(iv1_rfi_C, "HC1"))    # robust SE as in Stata 
iv1_rfi_Cr



############################
# Crime

# First Stage
iv1_fc_C <- lm(cristeros ~ bishop_dist100K +
                 Railways + telegraphs + 
                 rurales + morelos + mina + hidalgo + guerrero.x + french + 
                 dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                 pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                 loc_XVI + fran  + dom  + aug  + jes +  
                 arch_zone + triple_all + chichimeca , 
               data=bishops2)
summary(iv1_fc_C)
iv1_fc_Cr<-coeftest(iv1_fc_C, vcov = vcovHC(iv1_fc_C, "HC1"))    # robust SE as in Stata 
iv1_fc_Cr


# Second stage
iv1_sc_C<-ivreg(totalao2013 ~ cristeros + 
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca 
                | Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca + bishop_dist100K , 
                data=bishops2)  #, lag.instr=TRUE
summary(iv1_sc_C, vcov = sandwich, diagnostics = TRUE)
iv1_sc_Cr<-coeftest(iv1_sc_C, vcov = vcovHC(iv1_sc_C, "HC1"))    # robust SE as in Stata 
iv1_sc_Cr


# Reduced Form
iv1_rfc_C <- lm(totalao2013 ~ bishop_dist100K +
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca , 
                data=bishops2)
summary(iv1_rfc_C)
iv1_rfc_Cr<-coeftest(iv1_rfc_C, vcov = vcovHC(iv1_rfc_C, "HC1"))    # robust SE as in Stata 
iv1_rfc_Cr



############################
# DTOs

# First Stage
iv1_fd_C <- lm(cristeros ~ bishop_dist100K +
                 Railways + telegraphs + 
                 rurales + morelos + mina + hidalgo + guerrero.x + french + 
                 dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                 pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                 loc_XVI + fran  + dom  + aug  + jes +  
                 arch_zone + triple_all + chichimeca , 
               data=bishops2)
summary(iv1_fd_C)
iv1_fd_Cr<-coeftest(iv1_fd_C, vcov = vcovHC(iv1_fd_C, "HC1"))    # robust SE as in Stata 
iv1_fd_Cr


# Second stage
iv1_sd_C<-ivreg(total_all_dto10 ~ cristeros + 
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca 
                | Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca + bishop_dist100K , 
                data=bishops2)  #, lag.instr=TRUE
summary(iv1_sd_C, vcov = sandwich, diagnostics = TRUE)
iv1_sd_Cr<-coeftest(iv1_sd_C, vcov = vcovHC(iv1_sd_C, "HC1"))    # robust SE as in Stata 
iv1_sd_Cr


# Reduced Form
iv1_rfd_C <- lm(total_all_dto10 ~ bishop_dist100K +
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca , 
                data=bishops2)
summary(iv1_rfd_C)
iv1_rfd_Cr<-coeftest(iv1_rfd_C, vcov = vcovHC(iv1_rfd_C, "HC1"))    # robust SE as in Stata 
iv1_rfd_Cr




########################################################
# 6.8.2 Spatial 2SLS MODEL



si_C <-spreg(gini2 ~ cristeros + 
               Railways + telegraphs + 
               rurales + morelos + mina + hidalgo + guerrero.x + french + 
               dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
               pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
               loc_XVI + fran  + dom  + aug  + jes +  
               arch_zone + triple_all + chichimeca,
             data=bishops2, listw=nb_InvB,  
             initial.value = 0.2, model = "lag", het=TRUE)  
summary(si_C)

sink("../results/si_C.txt")
summary(si_C)
sink() 


sc_C<-spreg(totalao2013 ~ cristeros + 
              Railways + telegraphs + 
              rurales + morelos + mina + hidalgo + guerrero.x + french + 
              dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
              pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
              loc_XVI + fran  + dom  + aug  + jes +  
              arch_zone + triple_all + chichimeca,
            data=bishops2, listw=nb_InvB,  
            initial.value = 0.2, model = "lag", het=TRUE)  
summary(sc_C)

sink("../results/sc_C.txt")
summary(sc_C)
sink() 


sd_C<-spreg(total_all_dto10 ~ cristeros + 
              Railways + telegraphs + 
              rurales + morelos + mina + hidalgo + guerrero.x + french + 
              dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
              pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
              loc_XVI + fran  + dom  + aug  + jes +  
              arch_zone + triple_all + chichimeca,
            data=bishops2, listw=nb_InvB,  
            initial.value = 0.2, model = "lag", het=TRUE)  
summary(sd_C)

sink("../results/sd_C.txt")
summary(sd_C)
sink() 






########################################################
# 6.8.3 Spatial 2SLS MODEL



############################
# Inequality


# First stage
siv1_fhi_C_C<-spreg(cristeros ~ bishop_dist100K +
                      Railways + telegraphs + 
                      rurales + morelos + mina + hidalgo + guerrero.x + french + 
                      dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                      pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                      loc_XVI + fran  + dom  + aug  + jes +  
                      arch_zone + triple_all + chichimeca , 
                    data=bishops2, listw=nb_InvB, endog = NULL, instruments = NULL, 
                    lag.instr = TRUE, initial.value = 0.2, model = "lag", het=TRUE)
summary(siv1_fhi_C_C)

sink("../results/siv1_fhi_C_C.txt")
summary(siv1_fhi_C_C)
sink() 



## Impute data for gini in two observations
ginimean <- mean(bishops2$gini, na.rm = TRUE)
bishops2$gini2<-bishops2$gini
bishops2$gini2[is.na(bishops2$gini2)] <- ginimean


# Second stage
siv1_shi_C<-spreg(gini2 ~  Railways + telegraphs + 
                    rurales + morelos + mina + hidalgo + guerrero.x + french + 
                    dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                    pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                    loc_XVI + fran  + dom  + aug  + jes +  
                    arch_zone + triple_all + chichimeca , 
                  data=bishops2, listw=nb_InvB, endog = ~ cristeros, instruments = ~ bishop_dist100K, 
                  initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_shi_C)

# export
sink("../results/siv1_shi_C.txt")
summary(siv1_shi_C)
sink() 


# Reduced form
siv1_rfhi_C<-spreg(gini2 ~ bishop_dist100K + 
                     Railways + telegraphs + 
                     rurales + morelos + mina + hidalgo + guerrero.x + french + 
                     dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                     pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                     loc_XVI + fran  + dom  + aug  + jes +  
                     arch_zone + triple_all + chichimeca , 
                   data=bishops2, listw=nb_InvB,  
                   initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_rfhi_C)

# export
sink("../results/siv1_rfhi_C.txt")
summary(siv1_rfhi_C)
sink() 



############################
# Crime

# First stage
siv1_fhc<-spreg(cristeros ~ bishop_dist100K +
                  Railways + telegraphs + 
                  rurales + morelos + mina + hidalgo + guerrero.x + french + 
                  dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                  pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                  loc_XVI + fran  + dom  + aug  + jes +  
                  arch_zone + triple_all + chichimeca , 
                data=bishops2, listw=nb_InvB, endog = NULL, instruments = NULL, 
                lag.instr = TRUE, initial.value = 0.2, model = "lag", het=TRUE)
summary(siv1_fhc)

sink("../results/siv1_fhc.txt")
summary(siv1_fhc)
sink() 


# Second stage
siv1_shc_C<-spreg(totalao2013 ~  Railways + telegraphs + 
                    rurales + morelos + mina + hidalgo + guerrero.x + french + 
                    dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                    pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                    loc_XVI + fran  + dom  + aug  + jes +  
                    arch_zone + triple_all + chichimeca , 
                  data=bishops2, listw=nb_InvB, endog = ~ cristeros, instruments = ~ bishop_dist100K, 
                  initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_shc_C)

# export
sink("../results/siv1_shc_C.txt")
summary(siv1_shc_C)
sink() 


# Reduced form
siv1_rfhc_C<-spreg(totalao2013 ~ bishop_dist100K + 
                     Railways + telegraphs + 
                     rurales + morelos + mina + hidalgo + guerrero.x + french + 
                     dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                     pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                     loc_XVI + fran  + dom  + aug  + jes +  
                     arch_zone + triple_all + chichimeca , 
                   data=bishops2, listw=nb_InvB,  
                   initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_rfhc_C)

# export
sink("../results/siv1_rfhc_C.txt")
summary(siv1_rfhc_C)
sink() 



############################
# DTOs

# First stage
siv1_fhd_C_C<-spreg(cristeros ~ bishop_dist100K +
                      Railways + telegraphs + 
                      rurales + morelos + mina + hidalgo + guerrero.x + french + 
                      dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                      pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                      loc_XVI + fran  + dom  + aug  + jes +  
                      arch_zone + triple_all + chichimeca , 
                    data=bishops2, listw=nb_InvB, endog = NULL, instruments = NULL, 
                    lag.instr = TRUE, initial.value = 0.2, model = "lag", het=TRUE)
summary(siv1_fhd_C_C)

sink("../results/siv1_fhd_C_C.txt")
summary(siv1_fhd_C_C)
sink() 


# Second stage
siv1_shd_C<-spreg(total_all_dto10 ~  Railways + telegraphs + 
                    rurales + morelos + mina + hidalgo + guerrero.x + french + 
                    dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                    pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                    loc_XVI + fran  + dom  + aug  + jes +  
                    arch_zone + triple_all + chichimeca , 
                  data=bishops2, listw=nb_InvB, endog = ~ cristeros, instruments = ~ bishop_dist100K, 
                  initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_shd_C)

# export
sink("../results/siv1_shd_C.txt")
summary(siv1_shd_C)
sink() 


# Reduced form
siv1_rfhd_C<-spreg(total_all_dto10 ~ bishop_dist100K + 
                     Railways + telegraphs + 
                     rurales + morelos + mina + hidalgo + guerrero.x + french + 
                     dist_state_capital_log+elevation_log + dist_rail_k + gulf_3 + North_3 + Pacific_3 + 
                     pop_den_1921 + rural_p_1921 + cathol_p_1921  + agri_area_1928 + gov_fed_10K_1292 +  police_10K_1928 +  analfa_1921 +  mem_per_union_1929 + 
                     loc_XVI + fran  + dom  + aug  + jes +  
                     arch_zone + triple_all + chichimeca , 
                   data=bishops2, listw=nb_InvB,  
                   initial.value = 0.2, model = "lag", het=TRUE)  #, lag.instr=TRUE
summary(siv1_rfhd_C)

# export
sink("../results/siv1_rfhd_C.txt")
summary(siv1_rfhd_C)
sink() 










###############################################################
# End of script
###############################################################