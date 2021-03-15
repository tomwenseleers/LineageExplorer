# Function to plot the heatmap.

covidheat <- function(data, x = "DATE", y = "AGEGROUP", fill,
                      labfill = "Fill",
                      scale = c("absolute","relative")){
  x <- ensym(x)
  y <- ensym(y)
  fill <- ensym(fill)
  scale <- match.arg(scale)
  
  # Check if date
  if(!inherits(data[[!!x]], "Date")){
    data[[!!x]] <- as.Date(data[[!!x]])
  }
  
  # construct scale
  fillscale <- if(scale == "absolute"){
    scale_fill_viridis_c(option = "B")
  } else {
    scale_fill_gradient2(low = "darkblue",
                         mid = "white",
                         high = "darkred",
                         midpoint = 0)
  }
  
  ggplot(data, aes(x = !!x, y = !!y, fill = !!fill))+
    geom_raster() +
    theme_minimal() +
    labs(x = "Date", y = "Age group", fill = labfill) +
    scale_x_date(date_labels = "%b %d")
  
}
