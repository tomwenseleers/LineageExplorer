library(pals)
lineage_cols = glasbey(n)
pal.bands(lineage_cols)
pal.volcano(lineage_cols)
pal.zcurve(lineage_cols)
lineage_cols = kelly(n)
pal.bands(lineage_cols)
pal.volcano(lineage_cols)
pal.zcurve(lineage_cols)
lineage_cols = stepped(n) # OKish?
pal.bands(lineage_cols)
pal.volcano(lineage_cols)
pal.zcurve(lineage_cols)
lineage_cols = stepped2(n)
pal.bands(lineage_cols)
pal.volcano(lineage_cols)
pal.zcurve(lineage_cols)
lineage_cols = stepped3(n)
pal.bands(lineage_cols)
pal.volcano(lineage_cols)
pal.zcurve(lineage_cols)
lineage_cols = tableau20(n) # OKish?
pal.bands(lineage_cols)
pal.volcano(lineage_cols)
pal.zcurve(lineage_cols)
lineage_cols = watlington(n) # OKish?
pal.bands(lineage_cols)
pal.volcano(lineage_cols)
pal.zcurve(lineage_cols)

set.seed(5)
# lineage_cols = hsv(h = seq(0, 1, length=n), s = runif(n, 0.5, 1), v = runif(n, 0.5, 1))[1:n]
lineage_cols = hsv(h = seq(0, 1, length=n), s = runif(n, 0.5, 1), v = c(0.6, 0.75))[1:n]
# lineage_cols = lineage_cols[sample(1:n, size=n, replace=F)]
pal.bands(lineage_cols)
pal.volcano(lineage_cols)
pal.zcurve(lineage_cols)
lineage_cols = colorRampPalette(c("yellow2", "red2", "blue3", "magenta"))(n)
pal.bands(lineage_cols)
pal.volcano(lineage_cols)
pal.zcurve(lineage_cols)

banded2 = function (cols) sapply(1:length(cols), function (i) case_when((i %% 2)==1 ~ cols[[i]], 
                                                                        T ~ muted(cols[[i]], l=50, c=120) ) )
banded3 = function (cols) sapply(1:length(cols), function (i) case_when((i %% 3)==1 ~ cols[[i]],
                                                                        (i %% 3)==0 ~ muted(cols[[i]], l=50, c=130),
                                                                        (i %% 3)==2 ~ muted(cols[[i]], l=35, c=130) ) )
rainbow=rev(c(rgb(0.47,0.11,0.53),rgb(0.27,0.18,0.73),rgb(0.25,0.39,0.81),rgb(0.30,0.57,0.75),rgb(0.39,0.67,0.60),rgb(0.51,0.73,0.44),rgb(0.67,0.74,0.32),rgb(0.81,0.71,0.26),rgb(0.89,0.60,0.22),rgb(0.89,0.39,0.18),rgb(0.86,0.13,0.13))) 
rainbow2=rev(c(rgb(0.47,0.11,0.53),rgb(0.27,0.18,0.73),rgb(0.25,0.39,0.81),rgb(0.30,0.57,0.75),rev(c(rgb(0.39,0.67,0.60),rgb(0.51,0.73,0.44),rgb(0.67,0.74,0.32),rgb(0.81,0.71,0.26),rgb(0.89,0.60,0.22),rgb(0.89,0.39,0.18),rgb(0.86,0.13,0.13))))) 
rainbow3=c("purple","blue","green3","yellow2","red","magenta2")
lineage_cols_plot = c("grey65", banded2(colorRampPalette(rainbow)(n-2)), "magenta")
pal.bands(lineage_cols_plot)
pal.volcano(lineage_cols_plot)
pal.zcurve(lineage_cols_plot)

lineage_cols_plot = case_when(
  levels_VARIANTS_plot=="Other" ~ "grey65",
  levels_VARIANTS_plot=="Beta" ~ "green4",
  levels_VARIANTS_plot=="Alpha" ~ "#0085FF",
  levels_VARIANTS_plot=="Delta" ~ "mediumorchid",
  levels_VARIANTS_plot=="Omicron (BA.1)" ~ "red",
  levels_VARIANTS_plot=="Omicron (BA.2)" ~ "red3",
  levels_VARIANTS_plot=="Omicron (BA.2.12.1)" ~ "red4",
  levels_VARIANTS_plot=="Omicron (BA.4)" ~ "darkorange",
  levels_VARIANTS_plot=="Omicron (BA.4.6)" ~ "orange",
  levels_VARIANTS_plot=="Omicron (BA.5)" ~ "goldenrod",
  levels_VARIANTS_plot=="Omicron (BA.5.2)" ~ "yellow1",
  levels_VARIANTS_plot=="Omicron (BA.2.76)" ~ "magenta4",
  levels_VARIANTS_plot=="Omicron (BA.2.75)" ~ "magenta",
)

library(hues)
iwanthue(13, 0, 360, 54, 180, 27, 67, plot=TRUE)

lineage_cols_plot = c("grey65", polychrome(n)[-c(1:2)], "magenta")
pal.bands(lineage_cols_plot)
pal.volcano(lineage_cols_plot)
pal.zcurve(lineage_cols_plot)

lineage_cols_plot = case_when(
  levels_VARIANTS_plot=="Other" ~ "grey65",
  levels_VARIANTS_plot=="Beta" ~ "green4",
  levels_VARIANTS_plot=="Alpha" ~ "#0085FF",
  levels_VARIANTS_plot=="Delta" ~ "mediumorchid",
  levels_VARIANTS_plot=="Omicron (BA.1)" ~ "red",
  levels_VARIANTS_plot=="Omicron (BA.2)" ~ "red3",
  levels_VARIANTS_plot=="Omicron (BA.2.12.1)" ~ "red4",
  levels_VARIANTS_plot=="Omicron (BA.4)" ~ "darkorange",
  levels_VARIANTS_plot=="Omicron (BA.4.6)" ~ "orange",
  levels_VARIANTS_plot=="Omicron (BA.5)" ~ "goldenrod",
  levels_VARIANTS_plot=="Omicron (BA.5.2)" ~ "yellow1",
  levels_VARIANTS_plot=="Omicron (BA.2.76)" ~ "magenta4",
  levels_VARIANTS_plot=="Omicron (BA.2.75)" ~ "magenta",
)

lineage_cols = lineage_cols_plot[match(levels_VARIANTS,levels_VARIANTS_plot)]
