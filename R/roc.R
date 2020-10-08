#' Lehmann model for ROC curves.
#'
#' Relevant statistics for Lehmann models of ROC curves.
#'
#' @param formula A formula containing the covariates to be analyzed.
#' @param data A data.frame containing the variables named in the formula.
#'
#' @return The output is a list of relavant statistics. \code{theta},
#'  \code{var_theta}, \code{auc}, and \code{var_auc} are all numerics.
#'  \code{partial_auc} and \code{var_partial_auc} are both functions of the
#'  false positive rate. \code{roc} is a list containing numeric vectors
#'  \code{TPR} and \code{FPR}. \code{var_roc} is a function of the false
#'  positive rate. If the input formula contains concomitant covariates,
#'  the return list will contain functions. For example, \code{theta}
#'  would now be a function of the covariate values. The \code{apply_covariates}
#'  function should be used to input values of covariates into appropriate
#'  lehmann_roc objects.
#'
#' @examples
#' # Import a dataset
#' ovarian <- survival::ovarian
#' head(ovarian)
#'
#' # Create a Lehmann ROC model using only disease status.
#' l <- lehmann_roc(futime~resid.ds, ovarian)
#' summary(l)
#' plot(l)
#' l$theta
#' l$auc
#' l$partial_auc(0.5)
#'
#' # Create a Lehmann ROC model with age as a concomitant covariate.
#' l <- lehmann_roc(futime~resid.ds*age, ovarian)
#'
#' # The output is now a list of functions
#' summary(l)
#'
#' \dontrun{
#' l # errors
#' plot(l) # errors
#' }
#'
#' # Use the apply_covariates function
#' l_objects <- apply_covariates(l, list(age=21))
#' l2 <- l_objects[["21"]]
#' l2
#' l2$partial_auc(0.5)
#' plot(l2)
#'
lehmann_roc <- function(formula, data){
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
  }

  value <- list(theta=theta, var_theta=var_theta, auc=auc, var_auc=var_auc,
                partial_auc=partial_auc, var_partial_auc=var_partial_auc,
                roc=roc, var_roc=var_roc)
  attr(value, "class") <- "lehmann_roc"
  return(value)
}

#' Returns a list of lehmann_roc objects.
#'
#' This function is used when a user wants apply covariate values
#' to a lehmann_roc object. This function returns a list containing
#' n lehmann_roc objects, where n is the number of lists of covariates
#' provided.
#'
#' @param l A lehmann_roc object.
#' @param covariates A list of covariate values.
#' @param ... More lists of covariate values.
#'
#' @return A list containing lehmann_roc objects with the covariate
#' values applied. The labels of the list are the values of covariates,
#' separated by comams.
#'
#' @examples
#' \dontrun{
#' l <- # create lehmann_roc object
#' l_objects <- apply_covariates(l, list(cov_val_1=1, cov_val_2=2),
#' list(cov_val_1=2, cov_val_2=2))
#' l_1 <- l_objects[["1, 2"]]
#' l_2 <- l_objects[["2, 2"]]
#' l_1$theta
#' l_1$partial_auc(0.5)
#' plot(l_1)
#'
#' l_2$theta
#' l_2$partial_auc(0.5)
#' plot(l_2)
#' }
#'
apply_covariates <- function(l, covariates, ...){
  covariates <- list(covariates, ...)
  return_val <- list()
  for(cov_vals in covariates) {
    theta <- l$theta(cov_vals)
    var_theta <- l$var_theta(cov_vals)
    auc <- l$auc(theta)
    var_auc <- l$var_auc(theta, var_theta)
    partial_auc <- l$partial_auc(theta)
    var_partial_auc <- l$var_partial_auc(theta, var_theta)
    roc <- l$roc(theta)
    var_roc <- l$var_roc(theta, var_theta)

    lst <- list(theta=theta, var_theta=var_theta, auc=auc, var_auc=var_auc,
                partial_auc=partial_auc, var_partial_auc=var_partial_auc,
                roc=roc, var_roc=var_roc, covariates=cov_vals)
    attr(lst, "class") <- "lehmann_roc"
    return_val[[toString(cov_vals)]] <- lst
  }
  return(return_val)
}


# override default print function
print.lehmann_roc <- function(obj){
  cat("Theta:", obj$theta, "\n")
  cat("Variance of theta:", obj$var_theta, "\n")
  cat("AUC:", obj$auc, "\n")
  cat("Variance of AUC:", obj$var_auc, "\n")
}

# override default plot function
plot.lehmann_roc <- function(obj){
  plot(unlist(obj$roc[1]), unlist(obj$roc[2]),
       main="Lehmann ROC Curve", xlab="FPR", ylab="TPR")
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
    return(integrate(integrand, lower=0, upper=x0)$value)
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
