# Install and load the drc package
if(!require(drc)) {
    install.packages("drc")
    library(drc)
}

# Your data
concentration <- c(0.1, 0.37, 1.1, 3.3, 10, 30)
response <- c(1, 0.542936889, 0.344770451, 0.219157391, 0.145853623, 0.109774885)

# Fit a 4-parameter logistic model
model <- drm(response ~ concentration, fct = LL.4())

# Calculate IC50
IC50 <- ED(model, 50, interval = "delta")

# Output the IC50 value
print(IC50)
