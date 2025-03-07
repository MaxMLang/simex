#' Mixed SIMEX Method for Measurement Error and Misclassification
#'
#' This function implements the SIMEX method for models with both additive measurement error
#' and misclassification. It handles situations where one variable is misclassified and
#' another variable has additive measurement error.
#'
#' @author Max M. Lang
#' 
#' @aliases mixedsimex print.mixedsimex summary.mixedsimex print.summary.mixedsimex plot.mixedsimex predict.mixedsimex refit.mixedsimex
#'
#' @param model The naive model
#' @param SIMEXvariable Character vector specifying the names of variables with measurement error or misclassification
#' @param error.type Character vector indicating the type of error for each SIMEXvariable ("me" for measurement error, "mc" for misclassification)
#' @param measurement.error Measurement error variance-covariance matrix for variables with additive error
#' @param mc.matrix Misclassification probability matrix or list of matrices for variables with misclassification
#' @param lambda Vector of lambda values for the extrapolation step. Default is c(0.5, 1, 1.5, 2)
#' @param B Number of simulations for each lambda. Default is 100
#' @param fitting.method Method used for fitting the extrapolation. Default is "quadratic"
#' @param jackknife.estimation Specifying the extrapolation method for jackknife variance estimation. Can be set to FALSE if it should not be performed
#' @param asymptotic Logical, indicating if asymptotic variance estimation should be done
#' @param x Object of class 'mixedsimex'
#' @param digits Number of digits to be printed
#' @param object Object of class 'mixedsimex'
#' @param xlab Optional name for the X-Axis
#' @param ylab Vector containing the names for the Y-Axis
#' @param ask Logical. If TRUE, the user is asked for input, before a new figure is drawn
#' @param show Vector of logicals indicating for which variables a plot should be produced
#' @param newdata Optionally, a data frame in which to look for variables with which to predict
#' @param \dots Additional arguments passed to other functions
#'
#' @return An object of class "mixedsimex" containing the results of the SIMEX procedure.
#' \item{coefficients}{the corrected coefficients of the mixed SIMEX model,}
#' \item{SIMEX.estimates}{the estimates for every lambda,}
#' \item{model}{the naive model,}
#' \item{measurement.error}{the measurement error matrix,}
#' \item{mc.matrix}{the misclassification matrix,}
#' \item{B}{the number of iterations,}
#' \item{extrapolation}{the model object of the extrapolation step,}
#' \item{fitting.method}{the fitting method used in the extrapolation step,}
#' \item{residuals}{the residuals of the main model,}
#' \item{fitted.values}{the fitted values of the main model,}
#' \item{call}{the function call,}
#' \item{variance.jackknife}{the jackknife variance estimate,}
#' \item{extrapolation.variance}{the model object of the variance extrapolation,}
#' \item{variance.jackknife.lambda}{the data set for the extrapolation,}
#' \item{variance.asymptotic}{the asymptotic variance estimates,}
#' \item{theta}{the estimates for every B and lambda,}
#' ...
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with both measurement error and misclassification
#' # Generate data
#' set.seed(123)
#' n <- 200
#' x_true <- rnorm(n)
#' z_true <- rbinom(n, 1, 0.5)
#' y <- 1 + 2*x_true - 1.5*z_true + rnorm(n, 0, 0.5)
#' 
#' # Add measurement error to x
#' me.var <- 0.5
#' x_obs <- x_true + rnorm(n, 0, sqrt(me.var))
#' 
#' # Add misclassification to z
#' mc.matrix <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2)
#' z_obs <- z_true
#' for(i in 1:n) {
#'   if(z_true[i] == 0) {
#'     z_obs[i] <- sample(c(0,1), 1, prob = c(mc.matrix[1,1], mc.matrix[1,2]))
#'   } else {
#'     z_obs[i] <- sample(c(0,1), 1, prob = c(mc.matrix[2,1], mc.matrix[2,2]))
#'   }
#' }
#' 
#' # Create data frame
#' data <- data.frame(y = y, x_obs = x_obs, z_obs = factor(z_obs))
#' 
#' # Fit naive model
#' naive_model <- lm(y ~ x_obs + z_obs, data = data, x = TRUE)
#' 
#' # Apply mixed SIMEX
#' result <- mixedsimex(naive_model, 
#'                     SIMEXvariable = c("x_obs", "z_obs"),
#'                     error.type = c("me", "mc"),
#'                     measurement.error = me.var, 
#'                     mc.matrix = mc.matrix)
#' summary(result)
#' plot(result)
#' }
mixedsimex <- function(model, 
                       SIMEXvariable, 
                       error.type,
                       measurement.error = NULL, 
                       mc.matrix = NULL, 
                       lambda = c(0.5, 1, 1.5, 2), 
                       B = 100,
                       fitting.method = "quadratic", 
                       jackknife.estimation = "quadratic", 
                       asymptotic = TRUE) {
  
  # Check inputs
  fitting.method <- substr(fitting.method, 1, 4)
  if (!any(fitting.method == c("quad", "line", "nonl"))) {
    warning("Fitting Method not implemented. Using: quadratic", call. = FALSE)
    fitting.method <- "quad"
  }
  
  if (jackknife.estimation != FALSE)
    jackknife.estimation <- substr(jackknife.estimation, 1, 4)
  if (!any(jackknife.estimation == c("quad", "line", "nonl", FALSE))) {
    warning("Fitting Method (jackknife) not implemented. Using: quadratic",
            call. = FALSE)
    jackknife.estimation <- "quad"
  }
  
  if (!is.character(SIMEXvariable))
    stop("SIMEXvariable must be character", call. = FALSE)
  
  if (length(SIMEXvariable) != length(error.type))
    stop("SIMEXvariable and error.type must have the same length", call. = FALSE)
  
  if (!all(error.type %in% c("me", "mc")))
    stop("error.type must be either 'me' or 'mc'", call. = FALSE)
  
  if (any(lambda <= 0)) {
    warning("lambda should not contain 0 or negative values. 0 or negative values will be ignored",
            call. = FALSE)
    lambda <- lambda[lambda >= 0]
  }
  
  # Check model requirements
  if (class(model)[1] == "polr" && !any(names(model) == "Hessian"))
    stop("The option Hessian must be enabled in the naive model", call. = FALSE)
  
  if (!any(names(model) == "x") && asymptotic && class(model)[1] != "polr")
    stop("The option x must be enabled in the naive model for asymptotic variance estimation",
         call. = FALSE)
  
  if (class(model)[1] == "coxph" && asymptotic)
    stop("Asymptotic estimation is not supported for coxph models", call. = FALSE)
  
  if (class(model)[1] == "coxph" && is.null(model$model))
    stop("The option model = TRUE must be enabled for coxph models", call. = FALSE)
  
  # Handle Surv objects in coxph models
  if (class(model)[1] == "coxph" && grep("Surv\\(", names(model$model)[1]) == 1){
    timeEventMatrix <- as.matrix(model$model[[1]])
    timeName <- sub("Surv\\(","",strsplit(names(model$model)[1], ", ")[[1]][1])
    eventName <- sub("\\)","",strsplit(names(model$model)[1], ", ")[[1]][2])
    colnames(timeEventMatrix) <- c(timeName, eventName)
    model$model <- cbind(model$model, timeEventMatrix)
  }
  
  # Check measurement error matrix
  me_vars <- SIMEXvariable[error.type == "me"]
  if (length(me_vars) > 0 && is.null(measurement.error)) {
    stop("measurement.error must be provided for variables with measurement error", call. = FALSE)
  }
  
  if (!is.null(measurement.error)) {
    measurement.error <- as.matrix(measurement.error)
    SIMEXdata <- model$model
    
    # Handle homoscedastic vs heteroscedastic measurement error
    if (NROW(measurement.error) != NROW(SIMEXdata) && NROW(measurement.error) == 1){
      measurement.error <- matrix(measurement.error, nrow = NROW(SIMEXdata), ncol = NCOL(measurement.error), byrow = TRUE)
    }
    
    if (NROW(measurement.error) != NROW(SIMEXdata) && NROW(measurement.error) != 1){
      stop("NROW(measurement.error) must be either 1 or must take the number of rows of the data used.",
           call. = FALSE)
    }
    
    # Check for negative or zero measurement error
    if (any(measurement.error < 0))
      stop("measurement.error contains negative values", call. = FALSE)
    
    any0 <- (apply(measurement.error, 2, all.equal, current = rep(0, times = NROW(measurement.error))) == TRUE)
    if (sum(any0) > 0)
      stop("measurement.error is constant 0 in column(s) ", which(any0), call. = FALSE)
  }
  
  # Check misclassification matrix
  mc_vars <- SIMEXvariable[error.type == "mc"]
  if (length(mc_vars) > 0 && is.null(mc.matrix)) {
    stop("mc.matrix must be provided for variables with misclassification", call. = FALSE)
  }
  
  if (!is.null(mc.matrix)) {
    # Convert single matrix to list if needed
    if (is.matrix(mc.matrix) && length(mc_vars) == 1) {
      mc.matrix <- list(mc.matrix)
      names(mc.matrix) <- mc_vars
    }
    
    # Check mc.matrix format
    if (is.list(mc.matrix)) {
      if (length(mc.matrix) != length(mc_vars))
        stop("mc.matrix and misclassification variables do not match", call. = FALSE)
      
      # Check that variables are factors
      for (var in mc_vars) {
        if (!is.factor(model$model[[var]]))
          stop(paste("Misclassification variable", var, "must be a factor"), call. = FALSE)
      }
    } else if (!is.character(mc.matrix)) {
      stop("mc.matrix must be a matrix, list of matrices, or a function name", call. = FALSE)
    }
  }
  
  # Asymptotic variance limitations
  if (asymptotic == TRUE && 
      ((class(model)[1] != "glm" & class(model)[1] != "lm") || 
       (length(me_vars) > 0 && dim(unique(measurement.error))[1] != 1))) {
    stop("Asymptotic is only implemented for naive models of class lm or glm with homoscedastic measurement error.")
  }
  
  # Store call
  cl <- match.call()
  
  # Initialize parameters
  ncoef <- length(model$coefficients)
  if (class(model)[1] == "polr")
    ncoef <- ncoef + length(model$zeta)
  
  ndes <- dim(model$model)[1]
  
  if (class(model)[1] == "polr")
    p.names <- c(names(coef(model)), names(model$zeta))
  else
    p.names <- names(coef(model))
  
  nlambda <- length(lambda)
  estimates <- matrix(data = NA, nlambda + 1, ncoef) # +1 because "0" will be added
  theta <- matrix(data = NA, B, ncoef)
  colnames(theta) <- p.names
  theta.all <- vector(mode = "list", nlambda)
  
  # Initialize variance estimation
  if (jackknife.estimation != FALSE) {
    var.exp <- list()
    var.exp[[1]] <- extract.covmat(model)
  }
  
  if (asymptotic) {
    psi <- matrix(rep(0, ndes * ncoef), ncol = ncoef, nrow = ndes)
    psi <- resid(model, type = "response") * model$x
    PSI <- psi
    am <- list()
    a <- list()
    xi <- model$x
    dh <- rep(1, ndes)
    if (class(model)[1] == "glm")
      dh <- model$family$mu.eta(model$linear.predictors)
    for (k in 1:ndes)
      a[[k]] <- dh[k] * xi[k, ] %*% t(xi[k, ])
    a.mat <- matrix(unlist(a), nrow = length(a), byrow = TRUE)
    ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
    am[[1]] <- -ab / ndes
    a <- list()
  }
  
  # Assign the naive estimator
  if (class(model)[1] == "polr")
    estimates[1, ] <- c(model$coefficients, model$zeta)
  else
    estimates[1, ] <- model$coefficients
  
  # The Simulation step
  # Outer loop doing the simulations for each lambda
  for (i in 1:length(lambda)) {
    if (jackknife.estimation != FALSE)
      variance.est <- matrix(0, ncol = ncoef, nrow = ncoef)
    
    if (asymptotic) {
      psi <- matrix(0, ncol = ncoef, nrow = ndes)
      a <- list()
      for(k in 1:ndes)
        a[[k]] <- matrix(0, nrow = ncoef, ncol = ncoef)
    }
    
    # Inner loop, doing the simulations
    for (j in 1:B) {
      SIMEXdata <- model$model
      
      # Handle measurement error variables
      if (length(me_vars) > 0) {
        epsilon <- matrix(rnorm(n = NROW(SIMEXdata) * length(me_vars)),
                          ncol = length(me_vars),
                          nrow = NROW(SIMEXdata))
        
        # Add additional error based on lambda
        for (v in 1:length(me_vars)) {
          var_idx <- which(SIMEXvariable == me_vars[v])
          SIMEXdata[, me_vars[v]] <- SIMEXdata[, me_vars[v]] + 
            (sqrt(lambda[i]) * epsilon[, v] * measurement.error[, v])
        }
      }
      
      # Handle misclassification variables
      if (length(mc_vars) > 0) {
        for (v in 1:length(mc_vars)) {
          var_name <- mc_vars[v]
          
          # Get current categories
          categories <- levels(SIMEXdata[, var_name])
          n_categories <- length(categories)
          
          # Get misclassification matrix for this variable
          if (is.character(mc.matrix)) {
            # Use function-based misclassification
            mc_var_data <- data.frame(SIMEXdata[, var_name])
            colnames(mc_var_data) <- var_name
            SIMEXdata[, var_name] <- eval(call(mc.matrix, SIMEXdata, lambda[i]))[, var_name]
          } else {
            # Use matrix-based misclassification
            current_mc_matrix <- mc.matrix[[var_name]]
            
            # Check matrix dimensions
            if (nrow(current_mc_matrix) != n_categories || ncol(current_mc_matrix) != n_categories) {
              stop(paste("Misclassification matrix dimensions for", var_name, 
                         "do not match the number of categories"), call. = FALSE)
            }
            
            # Apply misclassification based on lambda
            if (lambda[i] > 0) {
              # Create lambda-dependent misclassification matrix
              mc_matrix_lambda <- current_mc_matrix
              
              # For lambda > 1, apply misclassification multiple times
              if (lambda[i] >= 1) {
                n_times <- floor(lambda[i])
                
                # Apply full misclassification n_times
                for (t in 1:n_times) {
                  mc_matrix_lambda <- mc_matrix_lambda %*% current_mc_matrix
                }
                
                # Apply partial misclassification for fractional part
                frac_part <- lambda[i] - n_times
                if (frac_part > 0) {
                  mc_matrix_frac <- (1 - frac_part) * diag(n_categories) + 
                    frac_part * current_mc_matrix
                  mc_matrix_lambda <- mc_matrix_lambda %*% mc_matrix_frac
                }
              } else {
                # For 0 < lambda < 1, interpolate between identity and mc.matrix
                mc_matrix_lambda <- (1 - lambda[i]) * diag(n_categories) + 
                  lambda[i] * current_mc_matrix
              }
              
              # Apply misclassification to each observation
              for (obs in 1:nrow(SIMEXdata)) {
                current_category <- SIMEXdata[obs, var_name]
                current_idx <- which(categories == current_category)
                
                # Sample new category based on misclassification probabilities
                new_idx <- sample(1:n_categories, 1, prob = mc_matrix_lambda[current_idx, ])
                SIMEXdata[obs, var_name] <- categories[new_idx]
              }
            }
          }
        }
      }
      
      # Update the model and calculate the estimate
      model.SIMEX <- update(model, data = data.frame(SIMEXdata))
      
      if (class(model)[1] == "polr")
        theta[j, ] <- c(model.SIMEX$coefficients, model.SIMEX$zeta)
      else
        theta[j, ] <- model.SIMEX$coefficients
      
      if (jackknife.estimation != FALSE) {
        variance.est <- variance.est + extract.covmat(model.SIMEX)
      }
      
      if (asymptotic) {
        xi <- model.SIMEX$x
        psi <- psi + (resid(model.SIMEX, type = "response") * xi)
        dh <- rep(1, ndes)
        if (class(model)[1] == "glm")
          dh <- model$family$mu.eta(model.SIMEX$linear.predictors)
        for (k in 1:ndes)
          a[[k]] <- a[[k]] - dh[k] * xi[k, ] %*% t(xi[k, ])
      }
    }
    
    # Taking the mean of the estimate -> SIMEX estimate
    estimates[i + 1, ] <- colMeans(theta)
    theta.all[[i]] <- theta
    
    # Variance estimation via the Jackknife
    if (jackknife.estimation != FALSE) {
      variance.est <- variance.est / B
      s2 <- cov(theta)
      var.exp[[i + 1]] <- variance.est - s2
    }
    
    if (asymptotic) {
      xiB <- psi / B
      PSI <- cbind(PSI, xiB)
      a.mat <- matrix(unlist(a), nrow = length(a), byrow = TRUE)
      ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
      am[[i + 1]] <- ab / (B * ndes)
    }
  }
  
  # Extrapolation step
  SIMEX.estimate <- vector(mode = "numeric", length = ncoef)
  colnames(estimates) <- p.names
  lambda <- c(0, lambda)
  
  # Fitting the extrapolation function
  switch(fitting.method,
         "quad" = extrapolation <- lm(estimates ~ lambda + I(lambda^2)),
         "line" = extrapolation <- lm(estimates ~ lambda),
         "nonl" = extrapolation <- fit.nls(lambda, p.names, estimates)
  )
  
  # Predicting the SIMEX estimate
  if (fitting.method == "nonl") {
    for (i in 1:length(p.names))
      SIMEX.estimate[i] <- predict(extrapolation[[p.names[i]]],
                                   newdata = data.frame(lambda = -1))
  } else {
    SIMEX.estimate <- predict(extrapolation, newdata = data.frame(lambda = -1))
  }
  
  # Jackknife Estimation
  if (jackknife.estimation != FALSE) {
    variance.jackknife <- matrix(unlist(var.exp), ncol = ncoef^2, byrow = TRUE)
    switch(jackknife.estimation,
           "quad" = extrapolation.variance <-
             lm(variance.jackknife ~ lambda + I(lambda^2)),
           "line" = extrapolation.variance <- lm(variance.jackknife ~ lambda),
           "nonl" = extrapolation.variance <-
             fit.nls(lambda, 1:NCOL(variance.jackknife), variance.jackknife)
    )
    variance.jackknife2 <- vector("numeric", ncoef^2)
    switch(jackknife.estimation,
           "nonl"= for(i in 1:NCOL(variance.jackknife))
             variance.jackknife2[i] <- predict(extrapolation.variance[[i]],
                                               newdata = data.frame(lambda = -1)),
           "quad"= variance.jackknife2 <- predict(extrapolation.variance,
                                                  newdata = data.frame(lambda = -1)),
           "line"= variance.jackknife2 <- predict(extrapolation.variance,
                                                  newdata = data.frame(lambda = -1))
    )
    variance.jackknife <- rbind(variance.jackknife2, variance.jackknife)
    variance.jackknife.lambda <- cbind(c(-1,lambda), variance.jackknife)
    variance.jackknife <- matrix(variance.jackknife[1, ], nrow = ncoef,
                                 ncol = ncoef, byrow = TRUE)
    dimnames(variance.jackknife) <- list(p.names, p.names)
  }
  
  # Asymptotic estimation
  if (asymptotic) {
    c11 <- cov(PSI)
    a11 <- diag.block(am)
    a11.inv <- solve(a11)
    sigma <- a11.inv %*% c11 %*% t(a11.inv)
    s <- construct.s(ncoef, lambda, fitting.method, extrapolation)
    d.inv <- solve(s %*% t(s))
    sigma.gamma <- d.inv %*% s %*% sigma %*% t(s) %*% d.inv
    g <- list()
    switch(fitting.method,
           "quad" = g <- c(1, -1, 1),
           "line" = g <- c(1, -1),
           "nonl" = for(i in 1:ncoef)
             g[[i]] <- c(-1, -(coef(extrapolation[[i]])[3] - 1)^-1,
                         coef(extrapolation[[i]])[2] / (coef(extrapolation[[i]])[3] - 1)^2)
    )
    g <- diag.block(g, ncoef)
    variance.asymptotic <- (t(g) %*% sigma.gamma %*% g) / ndes
    dimnames(variance.asymptotic) <- list(p.names, p.names)
  }
  
  # Creating class "mixedsimex"
  theta <- matrix(unlist(theta.all), nrow = B)
  theta.all <- list()
  for (i in 1:ncoef)
    theta.all[[p.names[i]]] <-
    data.frame(theta[, seq(i, ncoef * nlambda, by = ncoef)])
  
  z <- cbind(lambda, estimates)
  z <- rbind(c(-1, SIMEX.estimate), z) # returning the estimated values
  colnames(z) <- c("lambda", p.names)
  
  if(class(model)[1] == "polr")
    coefs <- z[1, -1][1:length(coef(model))]
  else
    coefs <- z[1, -1]
  
  erg <- list(
    coefficients = coefs, # SIMEX corrected coefficients
    SIMEX.estimates = z, # all thetas as a matrix
    lambda = lambda, # vector for the values for lambda
    model = model, # the naive model
    measurement.error = measurement.error, # measurement error matrix
    mc.matrix = mc.matrix, # misclassification matrix
    B = B, # number of Simulations
    extrapolation = extrapolation, # model of the extrapolation
    fitting.method = fitting.method, # which fitting method was used
    SIMEXvariable = SIMEXvariable,
    error.type = error.type,
    theta = theta.all,
    call = cl
  )
  
  class(erg) <- "mixedsimex"
  
  if (class(model)[1] == "polr")
    erg$zeta <- z[1,-1][(length(coef(model))+1):length(z[1,-1])]
  
  type <- 'response'
  if (class(model)[1] == "polr")
    type <- 'probs'
  if (class(model)[1] == "coxph")
    type <- 'lp'
  
  fitted.values <- predict(erg, newdata = model$model[, -1, drop = FALSE],
                           type = type)
  erg$fitted.values <- fitted.values
  
  if (class(model)[1] == "polr" || class(model)[1] == "coxph") {
    erg$residuals <- NULL
  } else if (is.factor(model$model[, 1])) {
    erg$residuals <-
      as.numeric(levels(model$model[, 1]))[model$model[, 1]] - fitted.values
  } else {
    erg$residuals <- model$model[, 1] - fitted.values
  }
  
  if (jackknife.estimation != FALSE) {
    erg$extrapolation.variance <- extrapolation.variance
    erg$variance.jackknife <- variance.jackknife
    erg$variance.jackknife.lambda <- variance.jackknife.lambda
  }
  
  if (asymptotic) {
    erg$PSI <- PSI
    erg$c11 <- c11
    erg$a11 <- a11
    erg$sigma <- sigma
    erg$sigma.gamma <- sigma.gamma
    erg$g <- g
    erg$s <- s
    erg$variance.asymptotic <- variance.asymptotic
  }
  
  return(erg)
}

#' @describeIn mixedsimex Print method for mixedsimex objects
#' @export
print.mixedsimex <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nNaive model:\n", deparse(x$model$call), "\n", sep = "")
  
  cat("\nSIMEX-Variables: ")
  me_vars <- x$SIMEXvariable[x$error.type == "me"]
  mc_vars <- x$SIMEXvariable[x$error.type == "mc"]
  
  if (length(me_vars) > 0) {
    cat("\n  Measurement error: ")
    cat(me_vars, sep = ", ")
  }
  
  if (length(mc_vars) > 0) {
    cat("\n  Misclassification: ")
    cat(mc_vars, sep = ", ")
  }
  
  cat("\nNumber of Simulations: ", paste(x$B), "\n\n", sep = "")
  
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2,
                  quote = FALSE)
  } else cat("No coefficients\n")
  
  if (length(x$zeta)) {
    cat("Intercepts:\n")
    print.default(format(x$zeta, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  
  cat("\n")
  return(invisible(x))
}

#' @describeIn mixedsimex Summary method for mixedsimex objects
#' @export
summary.mixedsimex <- function(object, ...) {
  if (class(object$model)[1] == "polr")
    est <- c(coef(object), object$zeta)
  else
    est <- coef(object)
  
  p.names <- names(est)
  est.table <- list()
  
  if (is.null(resid(object)))
    n <- object$model$n
  else
    n <- length(resid(object))
  
  p <- length(p.names)
  rdf <- n - p
  
  if (any(names(object) == "variance.jackknife")) {
    se <- sqrt(diag(object$variance.jackknife))
    tval <- est/se
    pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
    est.table[["jackknife"]] <- cbind(est, se, tval, pval)
    dimnames(est.table[["jackknife"]]) <- list(p.names, c("Estimate",
                                                          "Std. Error", "t value", "Pr(>|t|)"))
  }
  
  if (any(names(object) == "variance.asymptotic")) {
    se <- sqrt(diag(object$variance.asymptotic))
    tval <- est/se
    pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
    est.table[["asymptotic"]] <- cbind(est, se, tval, pval)
    dimnames(est.table[["asymptotic"]]) <- list(p.names, c("Estimate",
                                                           "Std. Error", "t value", "Pr(>|t|)"))
  }
  
  ans <- list()
  class(ans) <- "summary.mixedsimex"
  ans$coefficients <- est.table
  ans$residuals <- resid(object)
  ans$call <- object$call
  ans$B <- object$B
  ans$naive.model <- object$model$call
  ans$SIMEXvariable <- object$SIMEXvariable
  ans$error.type <- object$error.type
  ans$measurement.error <- object$measurement.error
  ans$mc.matrix <- object$mc.matrix
  ans$lambda <- object$lambda
  ans$fitting.method <- object$fitting.method
  
  return(ans)
}

#' @describeIn mixedsimex Print method for summary.mixedsimex objects
#' @export
print.summary.mixedsimex <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Mixed SIMEX Analysis\n")
  cat("--------------------\n")
  
  cat("Call:\n")
  print(x$call)
  
  cat("\nNaive model: \n")
  print(x$naive.model)
  
  cat("\nSIMEX variables:\n")
  me_vars <- x$SIMEXvariable[x$error.type == "me"]
  mc_vars <- x$SIMEXvariable[x$error.type == "mc"]
  
  if (length(me_vars) > 0) {
    cat("  Measurement error variables: ", paste(me_vars, collapse = ", "), "\n")
  }
  
  if (length(mc_vars) > 0) {
    cat("  Misclassification variables: ", paste(mc_vars, collapse = ", "), "\n")
  }
  
  cat("\nNumber of iterations:", x$B, "\n")
  cat("Lambda values:", paste(x$lambda, collapse = ", "), "\n")
  cat("Fitting method:", x$fitting.method, "\n\n")
  
  if (length(x$coefficients) > 0) {
    cat("Coefficients:\n")
    if (any(names(x$coefficients) == "jackknife")) {
      cat("Jackknife variance estimation:\n")
      printCoefmat(x$coefficients$jackknife, digits = digits)
    }
    
    if (any(names(x$coefficients) == "asymptotic")) {
      cat("\nAsymptotic variance estimation:\n")
      printCoefmat(x$coefficients$asymptotic, digits = digits)
    }
  } else {
    cat("No coefficient estimates available\n")
  }
  
  return(invisible(x))
}

#' @describeIn mixedsimex Plot method for mixedsimex objects
#' @export
plot.mixedsimex <- function(x,
                           xlab = expression((1 + lambda)),
                           ylab = colnames(b)[-1],
                           ask = FALSE,
                           show = rep(TRUE, NCOL(b) - 1), ...) {
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(...)
  
  if (ask)
    par(ask = TRUE)
  
  p.names <- names(coef(x))
  b <- x$SIMEX.estimates
  a <- seq(-1, max(b[, 1]), by = 0.01)
  d <- matrix(data = NA, nrow = length(a), ncol = NCOL(b) - 1)
  
  switch(x$fitting.method,
         quad = d <- matrix(predict(x$extrapolation, newdata = data.frame(lambda = a)), 
                           nrow = length(a), ncol = NCOL(b) - 1),
         line = d <- matrix(predict(x$extrapolation, newdata = data.frame(lambda = a)), 
                           nrow = length(a), ncol = NCOL(b) - 1),
         nonl = for (i in 1:length(p.names)) 
           d[, i] <- predict(x$extrapolation[[p.names[i]]], newdata = data.frame(lambda = a)))
  
  for (i in 2:NCOL(b)) {
    if (show[i - 1]) {
      plot(b[, 1] + 1, b[, i], xlab = xlab, ylab = ylab[i - 1], type = "n")
      points(b[-1, 1] + 1, b[-1, i], pch = 19)
      points(b[1, 1] + 1, b[1, i])
      lines(a[a > 0] + 1, d[a > 0, (i - 1)])
      lines(a[a < 0] + 1, d[a < 0, (i - 1)], lty = 2)
    }
  }
  
  invisible(x)
}

#' @describeIn mixedsimex Predict method for mixedsimex objects
#' @export
predict.mixedsimex <- function(object, newdata, ...) {
  new.object <- object$model
  new.object$coefficients <- object$coefficients
  
  if (class(new.object)[1] == "polr") {
    new.object$zeta <- object$zeta
    new.object$fitted.values <- object$fitted.values
  }
  
  if (missing(newdata)) {
    predict(new.object, ...)
  } else {
    predict(new.object, newdata = data.frame(newdata), ...)
  }
}

#' Refit method for mixedsimex objects
#'
#' @param object A mixedsimex object
#' @param fitting.method Method used for fitting the model. Default is "quadratic".
#' @param jackknife.estimation Method used for jackknife estimation. Default is "quadratic".
#' @param asymptotic Logical indicating whether to include asymptotic variance. Default is TRUE.
#' @param ... Additional arguments passed to the fitting function.
#'
#' @return A refitted mixedsimex object
#' @export
refit.mixedsimex <- function(object, fitting.method = "quadratic", 
                            jackknife.estimation = "quadratic",
                            asymptotic = TRUE, ...) {
  .refit(object, fitting.method = fitting.method,
         jackknife.estimation = jackknife.estimation, 
         asymptotic = asymptotic,
         allowed.fitting = c("quad", "line", "nonl"), 
         allowed.jackknife = c("quad", "line", "nonl", FALSE), ...)
}
