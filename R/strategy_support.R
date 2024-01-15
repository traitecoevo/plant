
##' Set a suitable hyperparameter function for chosen physiological model
##' @title Hyperparameters for FF16 physiological model
##' @param type Any strategy name as a string, e.g.: \code{"FF16"}.
##' @param parameters A parameters object
##' @rdname Hyperparameter_functions
##' @export
# if you update this function (even syntactic changes) update the function update_smc_support in the scaffolder
make_hyperpar <- function(type) {
  switch(type,
         FF16=make_FF16_hyperpar,
         FF16w=make_FF16w_hyperpar,
         FF16r=make_FF16r_hyperpar,
         K93=make_K93_hyperpar,
         stop("Unknown type ", type))
}

##' @rdname Hyperparameter_functions
##' @export
param_hyperpar <- function(parameters) {
  type <- attr(parameters$strategy_default, "class")
  switch(type,
         FF16_Strategy=FF16_hyperpar,
         FF16w_Strategy=FF16w_hyperpar,
         FF16r_Strategy=FF16r_hyperpar,
         K93_Strategy=K93_hyperpar,
         stop("Unknown type ", type))
}



##' @rdname Hyperparameter_functions
##' @export
# if you update this function (even syntactic changes) update the function update_smc_support in the scaffolder
hyperpar <- function(type) {
  switch(type,
         FF16=FF16_hyperpar,
         FF16w=FF16w_hyperpar,
         FF16r=FF16r_hyperpar,
         K93=K93_hyperpar,
         stop("Unknown type ", type))
}

##' @rdname make_environment
##' @export
environment_type <- function(type) {
  switch(type,
         FF16=sprintf("FF16_Env"),
         FF16w=sprintf("FF16_Env"),
         FF16r=sprintf("FF16_Env"),
         K93=sprintf("K93_Env"),
         stop("Unknown type ", type))
}

##' Make environment objects for a strategy
##' @param type Any strategy name as a string, e.g.: \code{"FF16"}.
##' @param parameters a object
##' @param ... other arguments passed through
##' @rdname make_environment
##' @export
make_environment <- function(type = NULL, parameters = NULL, ...) {
  
  if(!is.null(parameters)) {
    type = extract_RcppR6_template_types(parameters, "Parameters")[[1]][1]
#    plant_log_debug(sprintf("Creating default %s environment", type))
  }
    
  switch(type,
         FF16=FF16_make_environment(...),
         FF16w=FF16w_make_environment(...),
         FF16r=FF16r_make_environment(...),
         K93=K93_make_environment(...),
         stop("Unknown type ", type))
}

node_schedule_default <- function(p) {
  cl <- class(p)[[1]]
  switch(cl,
         "Parameters<FF16,FF16_Env>"=`node_schedule_default__Parameters___FF16__FF16_Env`,
         "Parameters<FF16w,FF16_Env>"=`node_schedule_default__Parameters___FF16w__FF16_Env`,
         "Parameters<FF16r,FF16_Env>"=`node_schedule_default__Parameters___FF16r__FF16_Env`,
         "Parameters<K93,K93_Env>"=`node_schedule_default__Parameters___K93__K93_Env`,
         stop("Unknown type: ", cl))(p)
}

make_node_schedule <- function(p) {
  cl <- class(p)[[1]]
  switch(cl,
         "Parameters<FF16,FF16_Env>"=`make_node_schedule__Parameters___FF16__FF16_Env`,
         "Parameters<FF16w,FF16_Env>"=`make_node_schedule__Parameters___FF16w__FF16_Env`,
         "Parameters<FF16r,FF16_Env>"=`make_node_schedule__Parameters___FF16r__FF16_Env`,
         "Parameters<K93,K93_Env>"=`make_node_schedule__Parameters___K93__K93_Env`,
                  stop("Unknown type: ", cl))(p)
}
