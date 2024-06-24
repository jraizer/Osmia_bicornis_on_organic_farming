# Function backward_selection
## Required packages
if(!"MuMIn" %in% installed.packages()) install.packages("MuMIn")
library(MuMIn)

backward_selection <- function(model) {
  current_model <- model
  terms <- attr(terms(model), "term.labels")
  
  model_info <- data.frame(Model = character(), 
                           AIC = numeric(), 
                           AICc = numeric(),
                           BIC = numeric(), 
                           LogLik = numeric(), 
                           stringsAsFactors = FALSE)
  
  while (length(terms) > 0) {
    aic_current <- AIC(current_model)
    aicc_current <- MuMIn::AICc(current_model)
    bic_current <- BIC(current_model)
    loglik_current <- logLik(current_model)
    
    model_info <- rbind(model_info, 
                        data.frame(Model = paste(terms, collapse = " + "), 
                                   AIC = aic_current, 
                                   AICc = aicc_current,
                                   BIC = bic_current, 
                                   LogLik = loglik_current))
    
    aicc_values <- sapply(terms, function(term) {
      reduced_formula <- as.formula(paste(". ~ . -", term))
      reduced_model <- update(current_model, reduced_formula)
      AICc(reduced_model)
    })
    aicc_min <- min(aicc_values)
    
    if (aicc_min < aicc_current) {
      term_to_remove <- names(which.min(aicc_values))
      reduced_formula <- as.formula(paste(". ~ . -", term_to_remove))
      current_model <- update(current_model, reduced_formula)
      terms <- attr(terms(current_model), "term.labels")
    } else {
      break
    }
  }
  
  model_info <- rbind(model_info, 
                      data.frame(Model = "Final Model", 
                                 AIC = AIC(current_model), 
                                 AICc = MuMIn::AICc(current_model),
                                 BIC = BIC(current_model), 
                                 LogLik = logLik(current_model)))
  
  return(list(final_model = current_model, model_info = model_info))
}

