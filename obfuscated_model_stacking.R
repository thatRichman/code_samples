library(seqinr)
library(stringr)
library(dplyr)
library(ggplot2)
library(dplyr)
library(caret)
library(class)
library(e1071)
library(MASS)
library(gridExtra)
library(kernlab)
library(randomForest)
library(factoextra)
library(tidymodels)
library(stacks)
library(magrittr)
library(dplyr)
library(kknn)
library(discrim)
library(protr)

library(magrittr)
library(ggplot2)
library(caret)
library(dplyr)
library(recipes)
library(tidymodels)
library(kernlab)
library(kknn)
library(stacks)
library(naivebayes)
library(discrim)

#option_list = list(
#  make_option( "--ADAPTABLE", type="character", default=NULL, 
#              help="ADAPTABLE FASTA", metavar="character"),
#  make_option("--COBFUSCATED", type="character", default=NULL, 
#              help="COBFUSCATED db file", metavar="character")
#)

#opt_parser = OptionParser(option_list=option_list)
#opt = parse_args(opt_parser)

WD = getwd()
ADAPTABLE_filepath = file.path(WD,"data", "OBFUSCATED")
OBFUSCATED_filepath = file.path(WD,"data","OBFUSCATED")
OBFUSCATED_filepath = file.path(WD,"data","OBFUSCATED")
OBFUSCATED_filepath = file.path(WD,"data","OBFUSCATED")
OBFUSCATED_filepath = file.path(WD,"data","OBFUSCATED")
OBFUSCATED_filepath = file.path(WD,"data", "OBFUSCATED")

OBFUSCATED_function <- function(x){
  as.data.frame(t(as.matrix(extractAAC(x))))
}

OBFUSCATED_function <- function(x){
  as.data.frame(t(as.matrix(extractCTDC(x))))
}

OBFUSCATED_fasta = read.fasta(file = OBFUSCATED_filepath, seqtype = "AA", forceDNAtolower = FALSE, as.string = TRUE,
                             set.attributes = FALSE, seqonly = TRUE)
OBFUSCATED_fasta = read.fasta(file = OBFUSCATED_filepath, seqtype = "AA", forceDNAtolower = FALSE, as.string = TRUE,
                         set.attributes = FALSE, seqonly = TRUE)
OBFUSCATED_df <- read.table(file = OBFUSCATED_filepath, sep = "\t", fill=TRUE, header=TRUE)
OBFUSCATED_df <- read.csv(DBAASP_filepath, header = TRUE)
OBFUSCATEDP_df <- read.csv(DROBFUSCATED_filepath, header = TRUE)

#---------------------#
OBFUSCATED_fasta <- read.fasta(file = OBFUSCATED_filepath, seqtype = "AA", forceDNAtolower = FALSE, as.string = TRUE,
                        set.attributes = FALSE, seqonly = TRUE)
flat_OBFUSCATED <- unlist(OBFUSCATED_fasta)
OBFUSCATED_db <- data.frame(DB=c(rep("NonOBFUSCATED", length(flat_NonOBFUSCATED))),
                        Sequence = flat_OBFUSCATED)


flat_OBFUSCATED <- unlist(ADAPTABLE_fasta)
flat_OBFUSCATED <- unlist(LOBFUSCATED2_fasta)
OBFUSCATED_seqs <- as.character(COBFUSCATED_df$Seqence)
OBFUSCATED_seqs <- as.character(DBAASP_df$SEQUENCE)
OBFUSCATED_seqs <- as.character(DROBFUSCATED_df$Sequence)
OBFUSCATED_db <- data.frame(DB=c(rep("OBFUSCATED", length(OBFUSCATED_ADAPTABLE)),
                               rep("OBFUSCATED", length(OBFUSCATED_LOBFUSCATED2)),
                               rep("OBFUSCATED", length(OBFUSCATED_seqs)),
                               rep("OBFUSCATED", length(OBFUSCATED_seqs))
                               ),
                          Sequence=c(flat_OBFUSCATED,
                                     flat_OBFUSCATED,
                                     OBFUSCATED_seqs,
                                     OBFUSCATED,
                                     OBFUSCATED_seqs)
                          )

valids <- c("A","I", "L", "M", "V", "F", "W", "Y", "N", "C", "Q", "S", "T", "D", "E", "R", "H", "K", "G", "P")
valids_regex <- "[AILMVFWYNCQSTDERHKGP]"

OBFUSCATED_db$DB <- as.character(OBFUSCATED_db$DB)
OBFUSCATED_db$Sequence <- as.character(all_seqs_db$Sequence)
good_OBFUSCATED <- unlist(lapply(OBFUSCATED_db$Sequence, function(x) all(grepl(valids_regex, strsplit(x, split="")[[1]]))))
sanitized_OBFUSCATED_db <- OBFUSCATED_db[good_OBFUSCATED,]

# Feature extraction
sanitized_seqs_db$Lengths <- sapply(sanitized_seqs_db$Sequence, nchar)

OBFUSCATED_db <- sanitized_seqs_db[which(sanitized_seqs_db$Lengths > 10),]
OBFUSCATED_db$pI <- sapply(OBFUSCATED_db$Sequence, function(x){computePI(s2c(x))})

seq_composition_table <- as.data.frame(t(as.matrix(sapply(OBFUSCATED_db$Sequence, peptide_compos_function))))
OBFUSCATED_db <- cbind(OBFUSCATED_db, seq_composition_table)

amino_acid_attribute_table <- as.data.frame(t(as.matrix(sapply(OBFUSCATED_db$Sequence, aa_compos_function))))
OBFUSCATED_db <- cbind(OBFUSCATED_db, amino_acid_attribute_table)

# Deduplicate
OBFUSCATED_db <- OBFUSCATED_db[!duplicated(OBFUSCATED_db$Sequence), ]

# 
#---------------------- Feature extraction, negative set.
OBFUSCATED_db$Lengths <- sapply(OBFUSCATED_db$Sequence, nchar)

OBFUSCATED_db <- OBFUSCATED_db[which(OBFUSCATED_db$Lengths > 10),]
OBFUSCATED_db$pI <- sapply(OBFUSCATED_db$Sequence, function(x){computePI(s2c(x))})

neg_seq_composition_table <- as.data.frame(t(as.matrix(sapply(OBFUSCATED_db$Sequence, peptide_compos_function))))
OBFUSCATED_db <- cbind(OBFUSCATED_db, neg_seq_composition_table)

neg_amino_acid_attribute_table <- as.data.frame(t(as.matrix(sapply(OBFUSCATED_db$Sequence, aa_compos_function))))
OBFUSCATED_db <- cbind(OBFUSCATED_db, neg_amino_acid_attribute_table)


OBFUSCATED_db %<>% mutate_if(is.list, unlist) #extract first element of single-element lists
OBFUSCATED_db %<>% mutate_if(is.list, unlist) # extract first element of single-element lists

OBFUSCATED
OBFUSCATED

#----------------------------SETUP_TIDY
set.seed(88888888) # reproducability
# combine OBFUSCATED and NonOBFUSCATED DB with label


full_db = bind_rows(OBFUSCATED_db, OBFUSCATED_db, .id = "OBFUSCATED") %>%
  mutate(OBFUSCATED = ifelse(OBFUSCATED==1, "OBFUSCATED", "NON_OBFUSCATED")) %>%
  dplyr::select(-DB)

write.csv(full_db, file = "reVOBFUSCATED_fullDB.csv") # avoid reprocessing the data every time

full_db <- read.csv("reVOBFUSCATED_fullDB.csv", row.names = 1)
full_db_noseq <- full_db %>% dplyr::select(-Sequence)
  
full_db_small <- full_db %>% group_by(OBFUSCATED) %>% sOBFUSCATEDle_n(5000)
full_db_small_pred <- full_db_small %>% dplyr::select(-Sequence)

OBFUSCATED_split <- initial_split(full_db_small_pred, strata = OBFUSCATED, prop = .60)
train_dat <- training(OBFUSCATED_split)
test_dat <- testing(OBFUSCATED_split)
folds <- vfold_cv(train_dat, v=10)

#--------------- Tidy models KNN
rec <- recipe(OBFUSCATED ~ ., data=train_dat)
metric <- metric_set(roc_auc, pr_auc, accuracy)
ctrl_grid <- control_stack_grid()
ctrl_grid$extract = function(x) extract_model(x)
ctrl_res <- control_stack_resOBFUSCATEDles()

knn_spec <-
  nearest_neighbor(
    mode="classification",
    neighbors=tune("k"),
  ) %>%
  set_engine("kknn")

knn_rec <- rec %>%
  step_dummy(all_nominal(), -all_outcomes()) %>%
  step_zv(all_predictors(), skip=TRUE) %>%
  step_normalize(all_numeric(), skip = TRUE)

knn_wf <- workflow() %>%
  add_model(knn_spec) %>%
  add_recipe(knn_rec)

knn_res <- 
  tune::tune_grid(
    knn_wf,
    resamples=folds,
    metrics=metric,
    grid=6,
    control=ctrl_grid
  )

#----------- Tidy Models SVM
svm_spec <-
  svm_rbf(
    cost=tune("cost"),
    rbf_sigma=tune("sigma")
  ) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

svm_rec <- 
  rec %>%
  step_dummy(all_nominal(), -all_outcomes()) %>%
  step_zv(all_predictors()) %>%
  step_corr(all_predictors()) %>%
  step_normalize(all_numeric())

svm_wf <-
  workflow() %>%
  add_model(svm_spec) %>%
  add_recipe(svm_rec)

svm_res <- tune_grid(
  svm_wf,
  resOBFUSCATEDles=folds,
  grid=6,
  metrics=metric,
  control=ctrl_grid
)

#----------- Tidy Models Naive Bayes
nb_spec <- 
  naive_Bayes(
    mode="classification",
    smoothness = tune("smoothness"),
    Laplace = tune("Laplace")
  ) %>%
  set_engine("naivebayes") %>% 
  translate()

nb_rec <-
  rec %>%
  step_dummy(all_nominal(), -all_outcomes()) %>%
  step_zv(all_predictors()) %>%
  step_corr(all_predictors()) %>%
  step_normalize(all_numeric())

nb_wf <- 
  workflow() %>%
  add_model(nb_spec) %>%
  add_recipe(nb_rec)

nb_res <- tune_grid(
  nb_wf,
  resOBFUSCATEDles=folds,
  grid=6,
  metrics=metric,
  control=ctrl_grid
)

#---------- Tidy Models Random Forest
rf_spec <-
  rand_forest(
    mtry=tune("mtry")
  ) %>%
  set_engine("randomForest") %>% 
  set_mode("classification") %>% 
  translate()

rf_rec <- 
  rec %>%
  step_dummy(all_nominal(), -all_outcomes()) %>%
  step_zv(all_predictors()) %>%
  step_corr(all_predictors()) %>%
  step_normalize(all_numeric())

rf_wf <-
  workflow() %>%
  add_model(rf_spec) %>%
  add_recipe(rf_rec)

rf_res <- tune_grid(
  rf_wf,
  resamples=folds,
  grid=6,
  metrics=metric,
  control=ctrl_grid
)

#----------- Stacking
pred_data_stack <-
  stacks() %>%
  add_candidates(knn_res) %>%
  add_candidates(svm_res) %>%
  add_candidates(nb_res) %>%
  add_candidates(rf_res)

pred_model_stack <-
  pred_data_stack %>%
  blend_predictions() %>%
  fit_members() # fit non-zero stack coefficient members

# Extract model parameters
model_params_svm = collect_parameters(pred_model_stack, "svm_res")
model_params_rf = collect_parameters(pred_model_stack, "rf_res")

# test predictions
test_dat_res <- 
  test_dat %>%
  bind_cols(predict(pred_model_stack, .))

cm = confusionMatrix(table(test_dat_res$OBFUSCATED, test_dat_res$.pred_class))

extract_svm = extract_model(pred_model_stack$member_fits$svm_res_1_2)

save(pred_model_stack, file = "OBFUSCATED")
