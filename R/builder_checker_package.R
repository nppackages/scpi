remove.packages("scpi")
library(devtools)
library(testthat)


setwd("D:/Dropbox/Research/Cattaneo/packages/R/scpi")
# Install and check
devtools::build()
devtools::install(upgrade = "never")
devtools::check(manual = FALSE)

# Run testthat
devtools::test()

# Prepare documentation
devtools::document()
shell('R CMD Rd2pdf . --output=man/manual.pdf --force --no-preview')
