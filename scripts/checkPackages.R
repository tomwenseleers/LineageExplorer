# Checks and loads the necessary packages
#----------------------------------------

toload <- c("shiny",
            "dplyr",
            "tidyr",
            "shinydashboard",
            "ggplot2",
            "patchwork",
            "shinyWidgets")

# library gives a warning when package cannot be found. 
# character.only forces library to only accept character values.
# This is needed to apply library over a character vector.
# logical.return forces library to return FALSE if 
# installation didn't succeed.

id <- suppressWarnings(
  sapply(toload, 
         library, 
         logical.return = TRUE,
         character.only = TRUE,
         quietly = TRUE,
         warn.conflicts = FALSE)
)

# Test installation and return the packages that need to
# be installed.
if(!all(id)){
  errmessage <- paste("To run the app, the following packages should be installed first:",
                      paste(toload[!id], collapse = "\n"), 
                      sep = "\n")
  stop(errmessage, call. = FALSE)
}
