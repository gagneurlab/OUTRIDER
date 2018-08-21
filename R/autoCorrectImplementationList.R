autoEncoderImplList <- list(
    R = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectR(ods, q, theta, BPPARAM=BPPARAM, ...)
	},
    PEER = function(ods, q, theta, BPPARAM=bpparam(), ...){
		peer(ods)
	},
    PEER_residual = function(ods, q, theta, BPPARAM=bpparam(), ...){
		peer(ods)
	},
    robustR = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectRCooksIter2(ods, q, theta, BPPARAM=BPPARAM, ...)
	},
    robustRM1 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectRCooksIter(ods, q, theta, BPPARAM=BPPARAM, ...)
	},
    robustRTheta = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectRCooksIterTheta(ods, q, theta, BPPARAM=BPPARAM, ...)
	},
    python = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectPython(ods, ...)
	},
    cooksR = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectRCooksIter3(ods, q, theta, BPPARAM=BPPARAM, ...)
	},
    pca = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectPCA(ods, q)
	},
    Rob1E3Pval25L10It10 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		Rob1E3Pval25L10It10(ods, q, debug=FALSE, BPPARAM=BPPARAM, ...)
	},
    Rob1E3PvalThetaMix100L10It10 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		Rob1E3PvalThetaMix100L10It10(ods, q, debug=FALSE, BPPARAM=BPPARAM, ...)
	},
    debug = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectRCooksIter2Debug(ods, q, theta, debug=FALSE, BPPARAM=BPPARAM, ...)
	},
    robustTheta = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectRCooksIter2Debug(ods, q, robust='iterative',
                noFirst=TRUE, internIter=100, modelTheta=TRUE, debug=FALSE, 
                initTheta=200, BPPARAM=bpparam())
	},
    mask25 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectRCooksMaskDebug(ods, q=q, debug=FALSE, ...)
	},
    maskCooksMix100L5It40 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		maskCooksMix100L5It40(ods, q=q, debug=FALSE, ...)
	},
    maskCooks25L5It40 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		maskCooks25L5It40(ods, q=q, debug=FALSE, ...)
	},
    RobTheta200 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectRCooksIter2Debug(ods, q=q, robust='iterative', 
                modelTheta=TRUE, initTheta=200, internIter=100, debug=FALSE, ...)
	},
    RobNoFTheta200 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		autoCorrectRCooksIter2Debug(ods, q=q, robust='iterative', 
                noFirst=TRUE, internIter=100, modelTheta=TRUE, debug=FALSE,
                initTheta=200)
	},
    robThetaFade200L20It25 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		robThetaFade200L20It25(ods, q, BPPARAM=BPPARAM)
	},
    robMix25L5I40 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		robMix25L5I40(ods, q, BPPARAM=BPPARAM)
	},
    RobPvalThetaMix100L5It40 = function(ods, q, theta, BPPARAM=bpparam(), ...){
		RobPvalThetaMix100L5It40(ods, q, BPPARAM=BPPARAM)
	},
    RobPval200L5It40 = function(ods, q, theta, ...){
		RobPval200L5It40(ods, q, ...)
	},
    RobPval25L5It40 = function(ods, q, theta, ...){
		RobPval25L5It40(ods, q, ...)
    },
    edPca = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=TRUE, ...)
    },
    edRand = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=FALSE, ...)
    },
    edRandRob = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=FALSE, robust=TRUE, ...)
    },
    edRandRobNF = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=FALSE, robust=TRUE, noFirstRob=TRUE, ...)
    },
    edPcaRob = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=TRUE, robust=TRUE, ...)
    },
    edPcaRobNF = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=TRUE, robust=TRUE, noFirstRob=TRUE, ...)
    },
    edPcaRobTMax500 = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=TRUE, robust=TRUE, noFirstRob=TRUE, 
                maxTheta=500, ...)
    },
    edPcaRobNfTMax250 = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=TRUE, robust=TRUE, noFirstRob=TRUE, 
                      maxTheta=250, ...)
    },
    edPcaRobNfTMax200 = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=TRUE, robust=TRUE, noFirstRob=TRUE, 
                      maxTheta=200, ...)
    },
    edRandRobNfTMax200 = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=FALSE, robust=TRUE, noFirstRob=TRUE, 
                      maxTheta=200, ...)
    },
    edRandRobNfTMax200FT = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=FALSE, robust=TRUE, noFirstRob=TRUE, 
                      maxTheta=200, firstTheta=TRUE, ...)
    },
    edPcaRobNfTMax200FT = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=TRUE, robust=TRUE, noFirstRob=TRUE, 
                      maxTheta=200, firstTheta=TRUE, ...)
    },
    edRandRobNfTMax250FT = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=FALSE, robust=TRUE, noFirstRob=TRUE, 
                      maxTheta=250, firstTheta=TRUE, ...)
    },
    edPcaRobNfTMax250FT = function(ods, q, theta, ...){
        autoCorrectED(ods, q, theta, usePCA=TRUE, robust=TRUE, noFirstRob=TRUE, 
                      maxTheta=250, firstTheta=TRUE, ...)
    },
    newED = function(ods, q, theta, ...){
        fitAutoencoder(ods, q, robust=TRUE, thetaRange=c(1, 250), 
                convergence=1e-5, minMu=0.00, noRobustLast=TRUE,  ...)
    }
)
