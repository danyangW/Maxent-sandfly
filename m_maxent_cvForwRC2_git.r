library(tidyverse)
library(dismo)
library(usdm)
library(rmaxent)
library(blockCV)
library(sf)
library(rJava)
library(here)

##### 1. Define functions ----
### 1.1 The function for model fitting of all data ----
fitFun <- function (data,         # input data contains presence-absence and env vars
                    spName,       # per species
                    predictorNames,    # names of predictors
                    vars,
                    nVar){
  
  # subset data per species
  subdata <- data %>% 
    filter ((species == spName)|is.na(species)) %>%  # bg points has no species name
    dplyr::select (X3035, Y3035, year, preab, all_of(predictorNames)) 

  # empty lists to store results  
  maList <- list()
  
  # fit per env var
  for (v in 1:nrow(vars)) {
    
    ma <- maxent(x=subdata %>% dplyr::select(all_of(vars %>% slice(v) %>% paste() )), 
                 p=subdata %>%  dplyr::select(preab), 
                 args=c("-J", "-P",  # -J = Jackknife, -P = response curve
                        paste0("betamultiplier=", regularization_multiplier),
                        paste0("outputformat=", output_format),
                        "linear=true", "quadratic=true", "hinge=true", "product=true","threshold=false")  
    )
    maList [[v]] <- ma
  }
  
  # export
  d_file <- sp_file(spName, "d", nVar, regularization_multiplier, "csv")
  if (!file.exists(d_file)) write_csv(subdata, d_file)
  
  write_rds(maList, sp_file(spName, "maList", nVar, regularization_multiplier, "rds"))
  
}

### 1.2 The function for model fitting and cv-testing ----
cvFitFUN <- function (data,         # input data contains presence-absence and env vars
                      inputfile,
                      spName,       # per species
                      predictorNames,    # names of predictors 
                      vars,
                      nVar,
                      k){     # fold for cross-validation 
  
  # subset data per species
  subdata <- data %>% 
    filter ((species == spName)|is.na(species)) %>%  # bg points has no species name
    dplyr::select (X3035, Y3035, year, preab, all_of(predictorNames)) 
  
  pts <- subdata %>%
    st_as_sf(coords = c("X3035", "Y3035"), crs=3035)
  
  # check autocorrelation distance
#  cv_spatial_autocor(x = pts, 
#                     column = "preab")

  # make folds
  set.seed(2024)
  scv <- cv_cluster(x = pts,
                    column = "preab", # optional: counting number of train/test records
                    k = k)
#  cv_plot(cv = scv, x = pts)
  
  folddata <- subdata %>% mutate (fold=scv[[2]])

  # empty lists to store results  
  E <- tibble (inputfile=NA,spName=NA,nPredictor=NA,mVar=NA, fold=NA,auc=NA,tss=NA)

  for (i in which(scv$records$test_1!=0)){                                                    
    
    # fit per env var
    for (v in 1:nrow(vars)) {

      # input datasets per train or test
      preabtrain <- folddata %>% filter(fold!=i) %>% dplyr::select(preab, all_of(vars %>% slice(v) %>% paste() )) %>% 
        drop_na() %>% dplyr::select(preab)
      envtrain <- folddata %>% filter(fold!=i) %>% dplyr::select(all_of(vars %>% slice(v) %>% paste() )) %>% drop_na()
      subtest <- folddata %>% filter(fold==i)
      envtest <- subtest %>% dplyr::select(all_of(vars %>% slice(v) %>% paste() ))%>% drop_na() 
      envtest_p <- subtest %>% filter (preab==1) %>% dplyr::select(all_of(vars %>% slice(v) %>% paste() ))%>% drop_na()  
      envtest_a <- subtest %>% filter (preab==0) %>% dplyr::select(all_of(vars %>% slice(v) %>% paste() ))%>% drop_na() 

      # model fitting
      mTry <- try(maxent(x=envtrain, 
                  p=preabtrain, 
                  args=c("-J", "-P",  # -J = Jackknife, -P = response curve
                         paste0("betamultiplier=", regularization_multiplier),
                         paste0("outputformat=", output_format),
                         "linear=true", "quadratic=true", "hinge=true", "product=true","threshold=false")), 
                  silent=T)
      
      m <- if (inherits(mTry, "try-error")) NA else mTry

      # cross validation
      eTry <- try (evaluate (m, 
                             p=envtest_p, 
                             a=envtest_a),
                   silent=T)
      
      if (inherits(eTry, "try-error")) {
        
        e1 <- tibble (inputfile=inputfile,spName=spName,nPredictor=nVar,
                      mVar=vars %>% slice(v) %>% paste(collapse = ", "),
                      fold=i,auc=NA,tss=NA)
      } else {
        
        e <- evaluate(model=m,
                      p=envtest_p,
                      a=envtest_a)
        
        e <- evaluate(model=m,
                      p=envtest_p,
                      a=envtest_a,
                      tr=threshold(e,stat="spec_sens"))    
        
        # calculate auc, tss on test data
        auc <- e@auc
        tss <- e@TPR + e@TNR-1
        
        # combine
        e1 <- bind_cols(inputfile, spName, nVar, 
                        vars %>% slice(v) %>% paste(collapse = ", "),
                        i,auc,tss) %>% 
          as_tibble() 
        names(e1) <- c("inputfile","spName","nPredictor","mVar","fold","auc","tss")
      }
      E <- bind_rows(E,e1)
    }
  }
  
  E <- E %>% filter(!is.na(inputfile))
  # save 
  efile_s <- read_csv(f_eval_sp(spName), show_col_types = FALSE) %>% bind_rows(., E)
  write_csv(efile_s, f_eval_sp(spName))
  
  eFile <- read_csv(f_eval, show_col_types = FALSE) %>% bind_rows(., E)
  write_csv(eFile, f_eval)
}


### 1.3 The function for computing effect size ----
sizeFun <- function (inputfile,
                     spName,
                     predictorNames,
                     vars,
                     nVar,
                     regularization_multiplier){   

  # load fitted models and train data 
  maList <- readRDS(sp_file(spName, "maList", nVar, regularization_multiplier, "rds"))
  d <- read_csv(sp_file(spName, "d", nVar, regularization_multiplier, "csv")) %>% 
    filter (preab==1)

  # tibble to store variable importance value
  vi <- tibble (mVar=apply(vars, 1, function(row) paste(row, collapse = ", ")),
    mSize=0)

  # calculate varSize per model (instead of per var)
  for (v in 1:length(maList)){

    tryLambdas <- try (rmaxent::parse_lambdas(maList[[v]])[[1]], silent=T)
    
    if (inherits(tryLambdas,"try-error")) {vi_m <- NA} else {
      
    # select one model 
    l_m <- rmaxent::parse_lambdas(maList[[v]])[[1]] %>% 
      as_tibble() %>%
      # make product vars two vars   
      mutate (var1 = case_when(type!="product"~ var,
                               type=="product"~ str_sub(var,start=1,end=str_locate(var,",")[,"start"]-1))) %>% 
      mutate (var2 = case_when(type=="product"~str_sub(var,start=str_locate(var,",")[,"start"]+1,end=nchar(var))))  
      # add threshold value for threshold feature                        
#      mutate (threshold=case_when(
#        type=="threshold" ~ str_sub(feature,start=2,
#                                    end=str_locate(feature,"<")[,"start"]-1)%>%
#          as.numeric())) 
      
    # combine train data and lambda file per model
    data_long <- purrr::map_dfr(seq_len(nrow(l_m)),~d) 
    l_m_long <- purrr::map_dfr(seq_len(nrow(d)),~l_m)
    combo <-  bind_cols(l_m_long, data_long)
    
    # read value of the second var in product feature
    combo1 <- combo %>%
      rowwise() %>%
      mutate(var2val = ifelse(is.na(var2),NA, get(var2))) %>% 
      mutate(var1val = get(var1))
  
    # calculate effect size by lambda*...
    combo2 <- combo1 %>% 
      mutate (mSize_perF_perP = ifelse(type == "linear", lambda,
                              ifelse(type == "quadratic", lambda*2*var1val,
#                                     ifelse(type == "threshold", ifelse(var1val < threshold, 0, lambda),
                                            ifelse(type == "forward_hinge", ifelse(var1val < min, 0, lambda),
                                                   ifelse(type == "reverse_hinge", ifelse(var1val < max, lambda, 0),
                                                          ifelse(type == "product", lambda*2*sign(var1val*var2val)*sqrt(abs(var1val*var2val)), NA))))))

    # per point sum coefficient across features as can be plus and minus
    combo3 <- combo2 %>% 
      group_by (X3035, Y3035, year) %>% 
      reframe(mSize_perP=sum(mSize_perF_perP))
    
    # average effect size (absolute value) across observations
    vi_m <- sum(abs(combo3$mSize_perP), na.rm=T)/nrow(d)
    
    # store in vi tibble
    vi$mSize[which(vi$mVar==vars %>% slice(v) %>% paste(.,collapse=", "))] <- vi_m
  }}
  
  vi <- vi %>% 
    arrange(desc(mSize))

   # combine
  mSize <- bind_cols (inputfile=rep(inputfile,nrow(vi)),
                      spName=rep(spName,nrow(vi)),
                      nPredictor=rep(nVar,nrow(vi)),
                      vi) %>% 
    as_tibble() 
  
  # export
  sizeFile <- read_csv(f_size, show_col_types = FALSE) %>% bind_rows(., mSize)
  write_csv(sizeFile, f_size)
}


### 1.4 Function to compute VIF ----
vifFUN <- function (data,         # input data contains presence-absence and env vars
                    inputfile,
                    spName,            # per species
                    predictorNames,     # names of predictors 
                    vars,
                    nVar) {
  
  # env data per species
  envdata <- data %>% 
    filter (species==spName) %>% 
    dplyr::select (all_of(predictorNames)) 

  # tibble to store variable importance value
  V <- tibble (mVar=apply(vars, 1, function(row) paste(row, collapse = ", ")),
    mVIF=0)

  # compute VIF per model
  for (v in 1:nrow(vars)) {
    
    envdata_v <- envdata %>% dplyr::select(all_of(vars %>% slice(v) %>% paste()))
    
    if (length(envdata_v)>=2) {
      
      mVIF <- vif(as.data.frame(envdata_v)) %>% 
        mutate (mVIF=max(VIF)) %>% 
        dplyr::select (mVIF) %>% 
        unique() %>% 
        pull()

    # store in vi tibble
    V$mVIF[which(V$mVar==vars %>% slice(v) %>% paste(collapse = ", "))] <- mVIF 
    }
  }
    
    # combine
    V <- bind_cols (inputfile=rep(inputfile,nrow(V)),
                    spName=rep(spName,nrow(V)),
                    nPredictor=rep(nVar,nrow(V)),
                    V) %>% 
      as_tibble() 
    
    # export
    vFile <- read_csv(f_vif, show_col_types = FALSE)
    vFile <- bind_rows(vFile, V)
    write_csv(vFile, f_vif)
}


##### 2. Input ----
### 2.1 config ----
dir_input  <- here("data", "input")
dir_output <- here("output", "cv_forwRC") # cv_forwRC stands for spatial cross validation, forward selection of all random combinations. 
dir_fig    <- here("figures")

# global result files
f_vif  <- file.path(dir_output, "mVIF_forwRC.csv")
f_size <- file.path(dir_output, "mSize_forwRC.csv")
f_eval <- file.path(dir_output, "e_forwRC.csv")

# per-species helpers
sp_dir     <- function(sp) file.path(dir_output, sp)
f_eval_sp  <- function(sp) file.path(sp_dir(sp), paste0("e_", sp, ".csv"))
sp_file    <- function(sp, prefix, nVar, reg, ext)
  file.path(sp_dir(sp), sprintf("%s_%s_%dvar_reg%s.%s", prefix, sp, nVar, reg, ext))

dir.create(dir_output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dir_fig, "reg"), recursive = TRUE, showWarnings = FALSE)

# parameter
output_format <- "cloglog"
k=5

### 2.2 start files ----
if (!file.exists(f_vif))  write_csv(tibble(inputfile=NA, spName=NA, nPredictor=NA, mVar=NA, mVIF=NA), f_vif)
if (!file.exists(f_size)) write_csv(tibble(inputfile=NA, spName=NA, nPredictor=NA, mVar=NA, mSize=NA), f_size)
if (!file.exists(f_eval)) write_csv(tibble(inputfile=NA, spName=NA, nPredictor=NA, mVar=NA, fold=NA, auc=NA, tss=NA), f_eval)

### 2.3 input data ----
inputfile <- c("train") # <---
data <- read_csv(file.path(dir_input, paste0(inputfile, ".csv"))) %>% 
  as_tibble() 
colSums(is.na(data))

# species names
spName <- intersect(sp_config$spName, unique(na.omit(data$species)))

### 2.4 define parameters for functions ----
sp_config <- tribble(
  ~spName,          ~regMul,
  "Ph.papatasi",     9,
  "Ph.perniciosus",  9,
  "S.minuta",        9,
  "Ph.sergenti",     9,
  "Ph.tobbi",        9,
  "Ph.mascittii",    9,
  "Ph.ariasi",       9,
  "S.dentata",       9,
  "Ph.neglectus",    7,
  "Ph.perfiliewi",   7,
  "Ph.simici",       7,
  "Ph.alexandri",    5)

# predictors ----
predictorBioclim <- c(                                                                
  "bio1_std","bio10_std","bio11_std","bio12_std","bio13_std","bio14_std","bio15_std",
  "bio16_std","bio17_std","bio18_std","bio19_std","bio2_std","bio3_std","bio4_std",
  "bio5_std","bio6_std","bio7_std","bio8_std","bio9_std",
  # "cmi_max_std","cmi_mean_std","cmi_min_std","cmi_range_std",#"fcf_std",
  "gdd0_std","gdd10_std", "gdd5_std",
  #"gddlgd0_std","gddlgd10_std","gddlgd5_std","gdgfgd0_std","gdgfgd10_std","gdgfgd5_std",
  "fgd_std","gsl_std","gst_std","lgd_std",  #"swb_year_std",         # <--- missing 10%-30% values
  "gsp_std",#"hurs_max_std","hurs_mean_std","hurs_min_std","hurs_range_std",
  "ngd0_std","ngd10_std","ngd5_std","npp_std",
#  "pet_penman_max_std","pet_penman_mean_std","pet_penman_min_std","pet_penman_range_std",
  "scd_std"#,"swb_std","vpd_max_std","vpd_mean_std","vpd_min_std","vpd_range_std"
  )
# catVar <- c("kg0","kg1","kg2","kg3","kg4","kg5")
# Maxent requires a minimum of three data points per category for categorical var

predictorLUH = c( "c3ann_std","c3per_std","c3nfx_std","c4ann_std","c4per_std","pastr_std","urban_std",
  "primf_std","primn_std","range_std","secdf_std","secdn_std","secma_std","secmb_std")

predictorBioclimM = c(
#  "cmi_01_std","cmi_02_std","cmi_03_std","cmi_04_std","cmi_05_std","cmi_06_std",
#  "cmi_07_std","cmi_08_std","cmi_09_std","cmi_10_std","cmi_11_std","cmi_12_std",
#  "hurs_01_std","hurs_02_std","hurs_03_std","hurs_04_std","hurs_05_std","hurs_06_std",
#  "hurs_07_std","hurs_08_std","hurs_09_std","hurs_10_std","hurs_11_std","hurs_12_std",
#  "pet_penman_01_std","pet_penman_02_std","pet_penman_03_std","pet_penman_04_std",
#  "pet_penman_05_std","pet_penman_06_std","pet_penman_07_std","pet_penman_08_std",
#  "pet_penman_09_std","pet_penman_10_std","pet_penman_11_std","pet_penman_12_std",
  "pr_01_std","pr_02_std","pr_03_std","pr_04_std","pr_05_std","pr_06_std",
  "pr_07_std","pr_08_std","pr_09_std","pr_10_std","pr_11_std","pr_12_std",
  "tas_01_std","tas_02_std","tas_03_std","tas_04_std","tas_05_std","tas_06_std",
  "tas_07_std","tas_08_std","tas_09_std","tas_10_std","tas_11_std","tas_12_std",
  "tasmax_01_std","tasmax_02_std","tasmax_03_std","tasmax_04_std","tasmax_05_std","tasmax_06_std",
  "tasmax_07_std","tasmax_08_std","tasmax_09_std","tasmax_10_std","tasmax_11_std","tasmax_12_std",
  "tasmin_01_std","tasmin_02_std","tasmin_03_std","tasmin_04_std","tasmin_05_std","tasmin_06_std",
  "tasmin_07_std","tasmin_08_std","tasmin_09_std","tasmin_10_std","tasmin_11_std","tasmin_12_std"	
#  "vpd_01_std","vpd_02_std","vpd_03_std","vpd_04_std","vpd_05_std","vpd_06_std",
#  "vpd_07_std","vpd_08_std","vpd_09_std","vpd_10_std","vpd_11_std","vpd_12_std"
  )
predictorStatic <- c("lith_std"#,"tpin3km_std","aspect_std","slope_std","pop_std"
                     )

predictorBio <- c("vert_S_std","vert_deltaS_std","vert_deltaSS_std","vert_I_std",
                  "mam_S_std","mam_deltaS_std","mam_deltaSS_std",
                  "bird_S_std","fbird_S_std","nfbird_S_std",
                  "bird_deltaS_std","fbird_deltaS_std","nfbird_deltaS_std",
                  "bird_deltaSS_std","fbird_deltaSS_std","nfbird_deltaSS_std"
#                  "bird_dS0090_std","bird_dS1000_std","bird_dS1510_std",
#                  "fbird_dS0090_std","fbird_dS1000_std","fbird_dS1510_std",    # <--- missing 20% values
#                  "nfbird_dS0090_std","nfbird_dS1000_std","nfbird_dS1510_std",
#                  "bird_dS0090#_std","bird_dS1000#_std","bird_dS1510#_std", 
#                  "fbird_dS0090#_std","fbird_dS1000#_std","fbird_dS1510#_std",
#                  "nfbird_dS0090#_std","nfbird_dS1000#_std","nfbird_dS1510#_std"
)

#predictorMon <- c("ws_std","d2m_std","t2m_std","src_std","stl1_std",
#                    "st_std","tp_std","swvl1_std","swv_std","sp_std",
#                    "ws_3m_std","d2m_3m_std","t2m_3m_std","src_3m_std","stl1_3m_std",
#                    "st_3m_std","tp_3m_std","swvl1_3m_std","swv_3m_std","sp_3m_std",
#                    "cmi_std","cmi_3m_std","hurs_std","hurs_3m_std","pet_penman_std",
#                    "pet_penman_3m_std","tasmax_std","tasmax_3m_std","tasmin_std",
#                    "tasmin_3m_std","vpd_std","vpd_3m_std","swb_year_std")

predictorNames = c(predictorBioclim, predictorLUH, predictorBioclimM, predictorStatic,predictorBio)

colnames(data)

# 2.5 species folder ----
for (s in spName){
  dir.create(sp_dir(s), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(f_eval_sp(s)))
    write_csv(tibble(inputfile=NA, spName=NA, nPredictor=NA, mVar=NA, fold=NA, auc=NA, tss=NA), f_eval_sp(s))
}

##### 3. Model fitting and model selection ----
### 3.1 Try regularization multiplier and decide ----

nVar <- 1

vars <- combn(predictorNames, nVar) %>% 
  t() %>% 
  as_tibble()
  
for (s in spName){
  
  regularization_multiplier <- sp_config$regMul[sp_config$spName == s]
  
  fitFun(data=data,         # input data contains presence-absence and env vars
         spName=s,       # per species
         predictorNames=predictorNames,    # names of predictors
         vars=vars,
         nVar=nVar)
  
  # read response curve
  maList <- readRDS(sp_file(s, "maList", 1, regularization_multiplier, "rds"))
  
  pdf(file.path(dir_fig, "reg", paste0("response_", s, "_reg", regularization_multiplier, ".pdf")))
  par (mfrow = c(3,3))
  for (n in 1:length(maList)){dismo::response(maList[[n]])}
  par(mfrow=c(1,1))
  dev.off()
}
# check response curve and decide which to use per species


### 3.2 Check folds for cv  =----
for (s in spName) {
  # subset data per species
  subdata <- data %>% 
    filter ((species == s)|is.na(species)) %>%  # bg points has no species name
    dplyr::select (X3035, Y3035, year, preab, all_of(predictorNames)) 
  
  pts <- subdata %>%
    st_as_sf(coords = c("X3035", "Y3035"), crs=3035)
  
  print(s)
  # make folds
  set.seed(2024)
  scv <- cv_cluster(x = pts,
                    column = "preab", # optional: counting number of train/test records
                    k = k)
}

### 3.3 Model fitting, cv, compute effect size ----
nVar_min <- 1
nVar_max <- 1

# parallel computing
#cl <- makeCluster(12)  # Use how many cores
#clusterExport(cl, "cvFitFUN")      # Export the function to the cluster

for (s in spName) {   
  
  # renew per species 
  gc()
  
  regularization_multiplier <- sp_config$regMul[sp_config$spName == s]
  
  nVar <- nVar_min 
  
  while (nVar<=nVar_max){

    # prepare env data
    vars <- combn(predictorNames, nVar) %>% 
      t() %>% 
      as_tibble()
    
    # first check vif, only var combi with vif <=10 enter model fitting 
    vifFUN (data=data,
            inputfile=inputfile,
            spName=s,           
            predictorNames=predictorNames,    
            vars=vars,
            nVar=nVar)
    
    vFile <- read_csv(f_vif)
    
    # select only not-multicolinear combi
    vars <- vFile %>% 
      filter(spName == s, nPredictor == nVar, mVIF <= 10) %>% 
      pull(mVar) %>% 
      str_split(", ", simplify = TRUE) %>% 
      as_tibble(.name_repair = ~ paste0("V", seq_along(.x)))

    # model fitting and cv
#    parLapply(cl, # Run the function in parallel
              cvFitFUN (data=data,                 # input data contains presence-absence and env vars
                        inputfile=inputfile,
                        spName = s,                # per species
                        predictorNames=predictorNames,  # names of predictors 
                        vars=vars,
                        nVar=nVar,
                        k=k)                       # k fold cross-validation
#    )
    if (nVar!=1){fitFun(data=data,     # input data contains presence-absence and env vars
                       spName=s,       # per species
                       predictorNames=predictorNames,    # names of predictors
                       vars=vars,
                       nVar=nVar)}
      
    # computing effect size per model
    sizeFun (spName=s,
             inputfile=inputfile,
             predictorNames=predictorNames,
             vars=vars,
             nVar=nVar,
             regularization_multiplier=regularization_multiplier)

    # go to next loop 
    nVar <- nVar+1
  }
}

#stopCluster(cl) # Stop the cluster after completion

#### END ----