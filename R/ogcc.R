#' O-glycan Cancer Classifier (OGCC)
#'
#' \code{ogcc} performs a group of cancer class predictions based on the expression of O-glycan-forming glycosyltransferase genes.
#'
#' @param x a numeric data frame or matrix of expression data where columns are the features and rows are cases/samples. The set of features must contain the minimum set of features for the specified model. See `Details`.
#' @param d a string indicating the RNA-Seq measurement type; either 'RSEM', 'RPKM' or 'TPM'. 'RSEM' by default.
#' @param model a character string indicating the model to perform the prediction task. It must be one of the following: "normal_tumor", "normal_tumor_reduced", "normal_tumor_types", "normal_tumor_types_reduced", "types" and "types_reduced".
#' @param output a character value indicating the form of prediction output of the model. It must be either "raw" or "probs".
#' @return a character vector in case of output = "raw" or a dataframe of class probalibilities in case of output = "probs"
#'
#' @author Ayman Abuelela; ayman.elkhodiery@kaust.edu.sa
#'
#' @details
#' \code{ogcc} classifies cancer samples according to one of a set of predefined classification models using O-glycan-forming genes (OGFGs).
#'
#' The minimum set of OGFGs required in any model can be displayed using \code{\link{getMSF}} function. The required feature for a specific model may differ according to the feature selection step of that model.
#'
#' To display the output labels the model was trained on and can predict, use the function \code{\link{getLabels}}. These models are not generic models. They were developed and trained on specific set of cancer classes. Going beyond the classes specified for each model may result in misclassification problems.
#'
#' \code{ogcc} can perform the following classification tasks based on OGFGTs expression profile:
#'
#' \enumerate{
#'       \item model \strong{\code{normal_tumor}}: predicts whether a sample is normal or cancer based on the expression profile of the OGFGs.
#'       \item model \strong{\code{types}}: predicts the cancer type.
#' }
#'
#' @importFrom stats predict
#'
#' @seealso \code{\link{getMSF}} and \code{\link{getLabels}}
#' @export

ogcc <- function(x, model = 'types', d = 'RSEM', output = 'raw') {
    ## check the class of x
    if (!is.matrix(x) & !is.data.frame(x)) {
        stop("ERROR: x must be a matrix or data frame")
    }

    ## check the class of x values
    apply(x, 2, function(x) if (is.numeric(x)) {} else {stop("ERROR: x must have numeric values")})

    ## check that the model name is correct
    model_names <- c('types', 'normal_tumor')
    if (!(model %in% model_names)) {
        stop("ERROR: model name is not found. It must be one of the following: 'normal_tumor' and 'types'.")
    }

    ## check that the data set has the minimum set of features
    msf <- getMSF(model = model, d = d)
    if (length(intersect(msf, colnames(x))) != length(msf)) {
        stop("ERROR: x does NOT have the minimum set of required features")
    }

    ## check d
    d_ <- c('RSEM', 'RPKM', 'TPM')
    if (!(d %in% d_)) {
        stop("ERROR: d must be either 'RPKM', 'RSEM' or 'TPM'.")
    }

    ## get model object
    model_obj <- models[[tolower(d)]][["pred"]][[model]]

    ## determine the preprocessing object
    pp_obj <- models[[tolower(d)]][["pp"]][[model]]

    ## preprocessing
    ppx <- predict(pp_obj, x)

    ## prediction
    predictions <- predict(model_obj, ppx, type = output)

    ## return output
    return(predictions)
}


#' Get minimum set of features (getMSF)
#'
#' \code{getMSF} prints the minimum set of features required by a model from \code{ogcc}.
#'
#' @param model is a prediction model name from \code{\link{ogcc}}. See \code{model} argument.
#' @param d a string indicating the RNA-Seq measurement type; either 'RSEM' or 'RPKM'. 'RSEM' by default.
#'
#' @return a character vector of the required features for a working model.
#'
#' @examples
#' msf <- getMSF("normal_tumor")
#' print(msf)
#'
#' @author Ayman Abuelela; ayman.elkhodiery@kaust.edu.sa
#'
#' @seealso \code{\link{ogcc}} and \code{\link{getLabels}}
#' @export

getMSF <- function(model = 'types', d = 'RSEM') {
    ## check model name
    model_names <- c('types', 'normal_tumor')
    if (!(model %in% model_names)) {
        stop("ERROR: model name is not found. It must be one of the following: 'normal_tumor' and 'types'.")
    }
    ## check d
    d_ <- c('RSEM', 'RPKM', 'TPM')
    if (!(d %in% d_)) {
        stop("ERROR: d must be either 'RPKM', 'RSEM' or 'TPM'.")
    }
    model_obj <- models[[tolower(d)]][["pred"]][[model]]
    return(model_obj$coefnames)
}

#' Get the output labels (getLabels)
#'
#' \code{getLabels} print the levels of the output variable a model has had trained on.
#' @param model is a prediction model name from \code{\link{ogcc}}. See \code{model} argument.
#'
#' @return a character vector of the class labels a model was trained on and can predict.
#'
#' @examples
#' labels <- getLabels("normal_tumor")
#' print(labels)
#'
#' @author Ayman Abuelela; ayman.elkhodiery@kaust.edu.sa
#'
#' @seealso \code{\link{ogcc}} and \code{\link{getMSF}}
#' @export

getLabels <- function(model = 'types', d = 'RSEM') {
    ## check model name
    model_names <- c('types', 'normal_tumor')
    if (!(model %in% model_names)) {
        stop("ERROR: model name is not found. It must be one of the following: 'normal_tumor' and 'types'.")
    }
    ## check d
    d_ <- c('RSEM', 'RPKM', 'TPM')
    if (!(d %in% d_)) {
        stop("ERROR: d must be either 'RPKM', 'RSEM' or 'TPM'.")
    }
    model_obj <- models[[tolower(d)]][["pred"]][[model]]
    return(model_obj$finalModel$classes)
}

#' @importFrom utils data
.onLoad <- function(libname, pkgname) {
    ## load modelsh
    data(model, package = pkgname, envir = parent.env(environment()))
}
