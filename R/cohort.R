##' @importFrom Rcpp evalCpp
##' @importFrom R6 R6Class
##' @useDynLib plant
NULL

cohort_schedule_max_time_default <- function(p) {
  cl <- class(p)[[1]]
  switch(cl,
         "Parameters<FF16,FF16_Env>"=`cohort_schedule_max_time_default__Parameters___FF16__FF16_Env`,
         "Parameters<FF16r,FF16r_Env>"=`cohort_schedule_max_time_default__Parameters___FF16r__FF16r_Env`,
         "Parameters<K93,K93_Env>"=`cohort_schedule_max_time_default__Parameters___K93__K93_Env`,
         
         
         stop("Unknown type: ", cl))(p)
}

cohort_schedule_default <- function(p) {
  cl <- class(p)[[1]]
  switch(cl,
         "Parameters<FF16,FF16_Env>"=`cohort_schedule_default__Parameters___FF16__FF16_Env`,
         "Parameters<FF16r,FF16r_Env>"=`cohort_schedule_default__Parameters___FF16r__FF16r_Env`,
         "Parameters<K93,K93_Env>"=`cohort_schedule_default__Parameters___K93__K93_Env`,

         stop("Unknown type: ", cl))(p)
}

make_cohort_schedule <- function(p) {
  cl <- class(p)[[1]]
  switch(cl,
         "Parameters<FF16,FF16_Env>"=`make_cohort_schedule__Parameters___FF16__FF16_Env`,
         "Parameters<FF16r,FF16r_Env>"=`make_cohort_schedule__Parameters___FF16r__FF16r_Env`,
         "Parameters<K93,K93_Env>"=`make_cohort_schedule__Parameters___K93__K93_Env`,         stop("Unknown type: ", cl))(p)
}
