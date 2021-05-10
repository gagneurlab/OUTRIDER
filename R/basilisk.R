
outrider2_env <- BasiliskEnvironment("outrider2_env", pkgname="OUTRIDER",
        packages=c("python=3.6.12", "anndata==0.7.5", "statsmodels==0.11.1",
                    "pandas==1.1.5", "numpy==1.19.5", "scikit-learn==0.23.1"),
        pip=c("tensorflow==2.4.1", "tensorflow-probability==0.12.1",
                    "nipals==0.5.2", "py_outrider==0.1.0"))