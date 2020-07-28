# LEHMANN ROC CLASS
# useful call: coxph(formula = Surv(time, status) ~ age + sex + age:sex, data = lung)


#' Lehmann model for ROC curves.
#'
#' Relevant statistics for a Lehmann model of an ROC curve.
#'
#' @param formula A formula containing the covariates to be analyzed.
#' @param data A data.frame containing the variables named in the formula.
#' @param cov_vals An optional list of covariate values.
#'
#' @return The output is a list of relavant statistics. \code{theta},
#'  \code{var_theta}, \code{auc}, and \code{var_auc} are all numerics.
#'  \code{partial_auc} and \code{var_partial_auc} are both functions of the
#'  false positive rate. \code{roc} is a list containing numeric vectors
#'  \code{TPR} and \code{FPR}. \code{var_roc} is a function of the false
#'  positive rate. If one or more covariates are named in the formula, all
#'  numeric values in the return list become functions of the covariates.
#'  If a list of covariates, \code{cov_vals}, is provided, the covariate values
#'  will be passed into the functions and the output will again contain
#'  numerics.
#'
#' @examples
#' # Import the ovarian dataset
#' ovarian <- survival::ovarian
#' head(ovarian)
#'
#' # Create a Lehmann ROC model using only disease status.
#' l <- lehmann_roc(futime~resid.ds, ovarian)
#' summary(l)
#' plot(l)
#' l$theta
#' l$auc
#'
#' # Create a Lehmann ROC model with age as a concomitant covariate.
#' l <- lehmann_roc(futime~resid.ds*age, ovarian)
#'
#' # The output is now a list of functions
#' summary(l)
#'
#' theta <- l$theta
#' theta_at_60 <- theta(list(age=60))
#' theta_at_60
#' l$auc(theta_at_60)
#' roc <- l$roc(theta_at_60)
#' plot(roc$FPR, roc$TPR)
#'
#' # Run the same model but include a list of covariate values.
#' l <- lehmann_roc(futime~resid.ds*age, ovarian, list(age=60))
#' summary(l)
#' plot(l)
#' l$theta
#' l$auc
lehmann_roc <- function(formula, data, cov_vals=NULL){
  #turn first term of formula into a Surv object
  formula = as.formula(paste("survival::Surv(", formula[2], ") ~", formula[3]))

  ph <- survival::coxph(formula, data)
  betas <- coef(ph)
  cov <- vcov(ph)

  # discard unnecessary betas and rows/cols in covariance matrix
  main_var <- names(betas)[1]
  ls <- get_needed(betas, cov, main_var)
  betas = ls$betas
  cov = ls$cov

  no_cov <- FALSE
  if(length(betas) == 1){
    no_cov <- TRUE
  }

  # if there are no covariates, lehmann_roc should return
  # numbers, not functions
  if (no_cov){
    theta <- exp(betas[[main_var]])
    var_theta <- variance_theta_no_cov(betas[[main_var]], cov[1])

    auc <- auc(theta)
    var_auc <- variance_auc(theta, var_theta)

    #these are functions
    partial_auc <- partial_auc(theta)
    var_partial_auc <- variance_partial_auc(theta, var_theta)

    roc <- roc(theta)
    var_roc <- variance_roc(theta, var_theta)
  } else {
    # these are all functions
    theta <- get_theta(betas, main_var)
    var_theta <- variance_theta(betas, cov, main_var)
    auc <- auc
    var_auc <- variance_auc
    partial_auc <- partial_auc
    var_partial_auc <- variance_partial_auc
    roc <- roc
    var_roc <- variance_roc
    #if a list of covariates is provided, apply them to the functions
    if (length(cov_vals == 0)){
      theta <- theta(cov_vals)
      var_theta <- var_theta(cov_vals)
      auc <- auc(theta)
      var_auc <- var_auc(theta, var_theta)
      partial_auc <- partial_auc(theta)
      var_partial_auc <- var_partial_auc(theta, var_theta)
      roc <- roc(theta)
      var_roc <- var_roc(theta, var_theta)
    }
  }


  # RETURN LEHMANN_ROC OBJECT
  value <- list(theta=theta, var_theta=var_theta, auc=auc, var_auc=var_auc,
                partial_auc=partial_auc, var_partial_auc=var_partial_auc,
                roc=roc, var_roc=var_roc)
  attr(value, "class") <- "lehmann_roc"
  return(value)
}


# override default print function
print.lehmann_roc <- function(obj){
  stopifnot(class(obj) == "lehmann_roc") # make sure obj is a lehmann_roc object
  cat("Theta:", obj$theta, "\n")
  cat("Variance of theta:", obj$var_theta, "\n")
  cat("AUC:", obj$auc, "\n")
  cat("Variance of AUC:", obj$auc, "\n")
}

# override default plot function
plot.lehmann_roc <- function(obj){
  stopifnot(class(obj) == "lehmann_roc")
  plot(obj$roc$FPR, obj$roc$TPR, main="Lehmann ROC Curve") # maybe add auc in text box
}

# Returns list of betas and covariance matrix without unnecessary variables
get_needed <- function(betas, cov, main_var){
  beta_list = list()
  for (label in labels(betas)){
    if (grepl(main_var, label, fixed=TRUE)){
      beta_list = c(beta_list, betas[label])
    }
  }
  needed_names = labels(beta_list)
  new_cov = cov[needed_names, needed_names]
  return(list(betas=beta_list, cov=new_cov))
}


get_theta <- function(betas, main_var){
  theta <- function(covariates){
    stopifnot(typeof(covariates) == "list")
    total_exponent = betas[[main_var]]
    betas = betas[2:length(betas)]
    for (beta_name in names(betas)){
      total_cov = 1
      for (cov_name in names(covariates)){
        if (cov_name != main_var & grepl(cov_name, beta_name, fixed=TRUE)){
          total_cov = total_cov * covariates[[cov_name]]
        }
      }
      total_exponent = total_exponent + total_cov*betas[[beta_name]]
    }
    return(exp(total_exponent))
  }
  return(theta)
}


variance_theta_no_cov <- function(beta, beta_variance){
  return(exp(2*beta)*beta_variance)
}


variance_theta <- function(betas, cov, main_var){
  variance <- function(covariates){
    get_formula <- function(){
      formula_string <- "~x1"
      count <- 2
      betas <- betas[2:length(betas)]
      for (beta_name in names(betas)){
        total_cov <- 1
        for (cov_name in names(covariates)){
          if (cov_name != main_var & grepl(cov_name, beta_name, fixed=TRUE)){
            total_cov <- total_cov * covariates[[cov_name]]
          }
        }
        formula_string <- sprintf("%s + %f*x%i", formula_string, total_cov, count)
        count <- count + 1
      }
      formula <- as.formula(formula_string)
      formula <- update(formula, ~exp(.))
      return(formula)
    }

    formula = get_formula()
    se <- msm::deltamethod(formula, unlist(betas), cov)
    return(se^2)
  }
}


auc <- function(theta){
  return((theta+1)^(-1))
}


variance_auc <- function(theta, variance_theta){
  return((theta+1)^(-4)*variance_theta)
}


partial_auc <- function(theta){
  integrand <- function(x){
    return((theta+1)^(-1)*x^(theta+1))
  }

  auc_to <- function(x0){
    return(integrate(integrand, lower=0, upper=x0))
  }

  return(auc_to)
}


variance_partial_auc <- function(theta, variance_theta){
  variance_at <- function(x0){
    outside_brackets <- (x0^(theta+1)/(theta+1))^2
    first_term <- (x0^theta*log(x0))^2*variance_theta/(x0^(theta+1))^2
    second_term <- variance_theta/(theta+1)^2
    third_term <- -2*x0^(theta+1)*log(x0)*variance_theta/(x0^(theta+1)*(theta+1))
    inside_brackets = first_term + second_term + third_term
    return(outside_brackets*inside_brackets)
  }
  return(variance_at)
}


roc <- function(theta){
  FPR <- seq(0, 1, 0.01)
  TPR <- FPR^theta
  return(list(FPR=FPR, TPR=TPR))
}


variance_roc <- function(theta, variance_theta){
  var_at <- function(x0){
    return((x0^theta*log(x0))^2*variance_theta)
  }
  return(var_at)
}
