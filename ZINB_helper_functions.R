logit_fit <- function(model_name) {
  Pi = 1
  Alpha = NA
  Mu = 1
  result = list("Method" = "Logit" , "Mu" = Mu, "Alpha" = Alpha, "Pi" = Pi)
  return(result)
}

NB_fit <- function(data.df, equation, model_name){
  m1 <- glm.nb(equation, data = data.df)
  Pi = 0
  Alpha = 1. / m1$theta
  Mu = predict(m1, type = "response")[[1]]
  result = list("Method" = "NB" , "Mu" = Mu, "Alpha" = Alpha, "Pi" = Pi)
  return(result)
}

ZINB_fit <- function(data.df, equation, model_name){
  m1 <- zeroinfl(equation, data = data.df, dist = "negbin")
  Pi <- predict(m1, type = "zero")[[1]]
  Alpha = 1. / m1$theta
  Mu = predict(m1, type = "count")[[1]]
  result = list("Method" = "ZINB" , "Mu" = Mu, "Alpha" = Alpha, "Pi" = Pi)
  return(result)
}

Fit_model <- function(data.df, gene, type_id){
  model_name = paste0(type_id,"_", gene)
  sub_df <- data.df[, gene]
  if (length(sub_df) >= 10) {
    if ((sum(sub_df !=0) / length(sub_df)) < 0.1){
    #if (all(sub_df == 0)){
      result <- logit_fit(model_name)
    } else if (all(sub_df != 0)){
      eq = paste0(as.character(gene), " ~ 1")
      result <- NB_fit(data.df, equation = as.formula(eq), model_name)
    } else {
      eq = paste0(as.character(gene), " ~ 1")
      result <- try(ZINB_fit(data.df, equation = as.formula(eq), model_name))
      if("try-error" %in% class(result)) {
        print(c("Was not able to fit ZINB for gene:", gene, "Zero%", sum(sub_df == 0) / length(sub_df) ))
        result <- logit_fit(model_name)
      }
    }
    return(result)
  } else {
    result = list("Method" = "Not enough data" , "Mu" = NA, "Alpha" = NA, "Pi" = NA)
    return(result)
  }
}

Get_contamination_factor <- function(cell_name, first_type, second_type){
  Beta_list[[c]][Beta_list[[c]][,"Type1"] == first_type & Beta_list[[c]][,"Type2"] == second_type,][,"B"]
}

Get_contamination_counts <- function(beta, yg){
  contamin_count1 <- list()
  contamin_count2 <- list()
  for (k in 0:yg){
    contamin_count1 <- c(contamin_count1, round(k / beta))
    contamin_count2 <- c(contamin_count2, round(( yg - k) / (1 - beta)))
  }
  result <- list("contamin_count1" = contamin_count1 , "contamin_count2"= contamin_count2)
  return(result)
}



Compute_pyg <- function(gene, yg, list1, list2, cell_name, first_type, second_type, count_threshold) {
    beta <- Get_contamination_factor(cell_name, first_type, second_type)
    con <- Get_contamination_counts(beta, yg)
    con1 <- unlist(con$contamin_count1) 
    con2 <- unlist(con$contamin_count2) 
    if (max(con1) > 2 * count_threshold | max(con2) > 2 * count_threshold) {
      return(0)
    } else {
      term1 <- unlist(Compute_likelihood_efficient(con1, Method = list1[[gene,"Method"]], Pi = list1[[gene,"Pi"]], Mu = list1[[gene,"Mu"]], Alpha = list1[[gene,"Alpha"]]))
      term2 <- unlist(Compute_likelihood_efficient(con2, Method = list2[[gene,"Method"]], Pi = list2[[gene,"Pi"]], Mu = list2[[gene,"Mu"]], Alpha = list2[[gene,"Alpha"]]))
      return(sum(term1 * term2))
    }
  
}

Compute_pure_pyg<- function(gene, yg, list, cell_name, first_type){
  unlist(Compute_likelihood_efficient(yg, Method = list[[gene,"Method"]], Pi = list[[gene,"Pi"]], Mu = list[[gene,"Mu"]], Alpha = list[[gene,"Alpha"]]))
}

Compute_likelihood_efficient <- function(contam_list, Method, Mu, Pi, Alpha) {
  if(Method == "NB"){
    #print("I am here1")
    #print(contam_list)
    n <- 1/Alpha
    prob <- n /(n + Mu)
    L <- list()
    for (x in contam_list){
      L <- c(L, dnbinom(x, n, prob, log = FALSE))
    }
  } else if (Method == "ZINB") {
    #print("I am here2")
    #print(contam_list)
    n <- 1/Alpha
    prob <- n /(n + Mu)  
    L <- list()
    for (x in contam_list){
      if (x == 0) {
        LNB0 <- Pi + (1 - Pi) * (1 + Alpha * Mu) ^ (-1 / Alpha)
        L <- c(L, LNB0)
      } else {
        LNB <- dnbinom(x, n, prob, log = FALSE)
        LNB <- (1 - Pi) * LNB
        L <- c(L, LNB)
      }
    }
    
  } else if (Method == "Logit") {
    #print("I am here3")
    #print(contam_list)
    L <- list()
    for (x in contam_list) {
      if (x == 0){
        L <- c(L, 1)
      } else {
        L <- c(L, 0)
      }
    }
  }
  return(L)
}
