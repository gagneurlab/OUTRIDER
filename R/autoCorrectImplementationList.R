autoEncoderImplList <- list(
    peer = function(ods, q, theta, BPPARAM=bpparam(), ...){
		peer(ods)
	},
    pca = function(ods, q, theta, BPPARAM=bpparam(), ...){
        autoCorrectPCA(ods, q)
    },
    ae = function(ods, q, theta, BPPARAM=bpparam(), ...){
        autoCorrectR(ods, q=q, BPPARAM=BPPARAM, ...)
	},
    newED = function(ods, q, theta, ...){
        fitAutoencoder(ods, q, lasso=FALSE, correctTheta='none', ...)
    },
    NLas_TCSf = function(ods, q, theta, ...){
        fitAutoencoder(ods, q, lasso=FALSE, correctTheta='sf', ...)
    },
    YLas_TCNo = function(ods, q, theta, ...){
        fitAutoencoder(ods, q, lasso=TRUE, correctTheta='none', ...)
    },
    YLas_TCSf = function(ods, q, theta, ...){
        fitAutoencoder(ods, q, lasso=TRUE, correctTheta='sf', ...)
    },
    NLas_TCSfed_NRob_NCR_NTC_YLAS_OWL = function(ods, q, theta, ...){
        fitAutoencoder(ods, q, lasso=TRUE, correctTheta=FALSE, useOptim=FALSE, ...)
    },
    ed_NRob_NCR_NTC_YLASENC_OWL = function(ods, q, theta, ...){
        fitAutoencoder(ods, q, lasso=TRUE, correctTheta=FALSE, useOptim=FALSE, L1encoder=TRUE, ...)
    },
    OUTRIDER = function(ods, q, ...){
        fitAutoencoder(ods, q, lasso=TRUE, useOptim=FALSE, ...)
    },
    NOWL_TCNo = function(ods, q, ...){
        fitAutoencoder(ods, q, lasso=TRUE, useOptim=FALSE, ...)
    },
    YOWL_TCNo = function(ods, q, ...){
        fitAutoencoder(ods, q, lasso=TRUE, useOptim=TRUE, ...)
    },
    NLas_TCNo = function(ods, q, ...){
        fitAutoencoder(ods, q, lasso=FALSE, ...)
    }
)

names(autoEncoderImplList) <- tolower(names(autoEncoderImplList))
