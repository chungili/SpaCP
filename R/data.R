#' Training Dataset
#'
#' A dataset containing survival analysis data for model training.
#'
#' @format A list object containing the following components:
#' \describe{
#'   \item{Yij}{Survival time for each subject}
#'   \item{Xij}{Covariate matrix}
#'   \item{delta}{Event indicator (1 = event occurred, 0 = censored)}
#'   \item{Mi}{Spatial location of each subject}
#'   \item{d2}{Distance matrix between locations}
#' }
"train.dt"

#' Testing Dataset
#'
#' A dataset containing survival analysis data for model testing.
#'
#' @format the testing data
#' \describe{
#'   \item{Yij}{survial time}
#'   \item{Xij}{covaiates}
#'   \item{delta}{event}
#'   \item{Mi}{Location}
#'   \item{d2}{distance matrix}
#' }
"test.dt"