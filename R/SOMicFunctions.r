par.month.to.day = function(fit.par){
  fit.par[3:8]   = fit.par[3:8]   * (12/365)
  fit.par
}

par.day.to.month = function(fit.par){
    fit.par[3:8] = fit.par[3:8] / (12/365)
}

with.dt = function(dt, expr){
  # function to evaluate list of expressions on data.table, sequentially
  for (j in 1:length(expr)) set(dt, NULL, names(expr)[j], dt[, eval(expr[[j]])])
}

#' RothC temperature response function
#'
#' Calculates temperature rate modifier
#'
#' @param Temp Temperature in degrees C.
#' @return temperature rate modifier.
#' @author Dominic Woolf
#' @export
fT.RothC = function (Temp) {
  ifelse(Temp >= -18.3, 47.9/(1 + exp(106/(Temp + 18.3))), 0)
}

#' PrimC temperature response function
#'
#' Calculates temperature rate modifier
#'
#' @param Temp Temperature in degrees C.
#' @return temperature rate modifier.
#' @author Dominic Woolf
#' @export
fT.PrimC = function(temp, method = 'Century2', t.ref = 30){
  # function to calculate temperature dependence of decomposition rate
  switch(method,
    RothC = fT.RothC(temp),
    Century1 = fT.Century1(temp) * fT.RothC(t.ref) / fT.Century1(t.ref),
    Century2 = fT.Century2(temp) * fT.RothC(t.ref) / fT.Century2(t.ref),
    stop('Unrecognised temperature function'))
}


#' Estimate saturated water holding capacity
#'
#' Calculates saturation point
#'
#' @param Sand Sand in percent.
#' @param Silt Silt in percent.
#' @param Clay Clay in percent.
#' @return saturated water holding capacity (m/m).
#' @author Dominic Woolf
#' @export
wsat = function(sand, silt, clay) 0.6658*silt + 0.1567*sand - 0.0079*silt^2 - 12.31121/sand -
  6.4756*log(sand) - 0.0038*clay*silt + 0.0038*clay*sand - 0.0042*silt*sand + 52.7526

#' Estimate field water holding capacity
#'
#' Calculates field capacity
#'
#' @param Sand Sand in percent.
#' @param Silt Silt in percent.
#' @param Clay Clay in percent.
#' @return Field water holding capacity (m/m).
#' @author Dominic Woolf
#' @export
wfield = function(sand, silt, clay) 118.932*clay + 119.0866*silt + 119.1104*sand + 162.31731/clay -
  46.21921/silt-5.12991/sand  + 18.1733*log(clay) + 0.0013*clay*silt + 0.0022*silt*sand - 11939.3493

#' Estimate wilting point
#'
#' Calculates wilting point
#'
#' @param Sand Sand in percent.
#' @param Silt Silt in percent.
#' @param Clay Clay in percent.
#' @return Wilting point (m/m).
#' @author Dominic Woolf
#' @export
wwilt = function(sand, silt, clay) 92.3851 - 1.5722*silt - 0.5423*sand - 0.0072*clay^2 + 0.0072*silt^2 -
  0.0059*sand^2 + 160.14591/clay  +  6.60011/sand + 0.0022*clay*silt - 0.0039*clay*sand

#' Convert rasters to 0 to 360 longitudes
#'
#' Provides the inverse function of raster::rotate
#'
#' @param Sand Sand in percent.
#' @param Silt Silt in percent.
#' @param Clay Clay in percent.
#' @return Unrotated raster.
#' @author Dominic Woolf
#' @export
unrotate = function(x) { # inverse of rotate: convert to 0 to 360 longitudes
  raster::shift(rotate(raster::shift(x, 180)), 180)
}

#' Initialise Draw a pretty map
#'
#' Draws an SOC map
#'
#' @param data_path Path to data.
#' @param pattern filename pattern
#' @param full.names full names
#' @param recursive look in subdirectories too?
#'
#' @return A data.frame for the somic model.
#' @author Dominic Woolf
#' @note This is a very simple function.
#' @rdname get_data
#' @export
g.map = function (r){
  myPalette = colorRampPalette(rev(brewer.pal(11, "Spectral")))
  sf = scale_fill_gradientn(
    colours = myPalette(50),
    trans = "log",
    breaks = (c(1,2,5,10,20,40,80)),
    limits = (c(0.99,80)),
    na.value = NA,
    name = expression("SOC (kg m"^-2*")"))
  t = theme(axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position=c(0.42,-0.1),
    legend.title.align=1,
    legend.title = element_text(size=12),
    legend.text = element_text(size = 12),
    legend.direction="horizontal",
    legend.justification = c(0.4,0),
    plot.margin = unit(c(0,0,0,0), "cm"))
  
  r[r<1] = 1
  ggplot((r)) +
    geom_path(data = grat, mapping=aes(long,lat, group=group), size=0.3, linetype=2, color='grey50') +
    geom_raster(aes(fill=value)) +
    geom_path(data = wmap_df, mapping = aes(long, lat, group=group), size=0.3) +
    geom_path(data=bbox_df, mapping=aes(long,lat, group=group), size=0.3) +
    guides(fill= guide_colorbar(barwidth=12, title.vjust=1, nbin=500, draw.ulim = FALSE, draw.llim = FALSE)) +
    coord_equal() +
    map_cols +
    map_theme
}


#' Initialise Somic data from files
#'
#' Initialises model inputs from files
#'
#' @param data_path Path to data.
#' @param pattern filename pattern
#' @param full.names full names
#' @param recursive look in subdirectories too?
#'
#' @return A data.frame for the somic model.
#' @author Dominic Woolf
#' @note This is a very simple function.
#' @rdname get_data
#' @export
get_data = function(data_path) {
  files <- list.files(path=data_path, pattern="*.socdat", full.names=T, recursive=FALSE)
  experiments <<- basename(file_path_sans_ext(files))
  all.data <- setDT(ldply(files, function(fn)  data.frame(read.table(fn, header = T, sep = ","),exp=basename(file_path_sans_ext(fn)))))
  all.data$added.bio <- as.double(all.data$added.bio)
  all.data <- merge (all.data, exp.const, by="exp")
  setnames(all.data, c("init.soc", "month"), c("soc", "time"))
  all.data
}

#' Find the Mode (most likely value)
#'
#' Finds mode.
#'
#' @param x A vector.
#'
#' @return The mode.
#' \item{cdiff}{The mode.}
#' \item{x}{A vector.}
#' @author Dominic Woolf
#' @rdname cdiff
#' @export
Mode <- function(x, na.rm = FALSE) {
  if (na.rm) x = x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Initialise Somic data
#'
#' Initialises model inputs by adding required columns to input data.frame using default values
#'
#' @param soc.data A data frame that already has a time column.
#'
#' @return The data.table but with more stuff filled in.
#' \item{soc.data}{The data.frame that already has a time column..}
#' @author Dominic Woolf
#' @note This is a very simple function.
#' @rdname initialise.somic.data
#' @export
initialise.somic.data <- function(soc.data, init.soc = 0.0, init.soc.d13c = 0.0, init.cover=1L) {
  if (!is.data.table(soc.data)) setDT(soc.data)
  if (!('time' %in% names(soc.data))) stop('time column required')
  if (!('temp' %in% names(soc.data))) stop('temp column required')
  if (!('precip' %in% names(soc.data))) stop('precip column required')
  if (!('pet' %in% names(soc.data))) stop('pet column required')
  if (!('soc' %in% names(soc.data))) soc.data[, soc := init.soc]
  soc.data[, ipm := soc * 0.081]
  soc.data[, spm := soc * 0.001]
  soc.data[, doc := 0.0]
  soc.data[, mb := soc * 0.0331]
  soc.data[, mac := soc * (1-0.001-0.081-0.0331)]
  soc.data[, atsmd := 0.0]
  soc.data[, h2o := 0.0]
  soc.data[, a := fT.PrimC(temp, method = 'Century1', t.ref = 28)]
  if (!('cover' %in% names(soc.data))) soc.data[, cover := init.cover]
  soc.data[, c := ifelse (cover==0, 1, 0.6)]
  soc.data[, mic := 1]
  soc.data[, added.doc := 0.0]
  soc.data[, velocity := 0.0]
  soc.data[, add_14c_age := 0.0]
  if (!('soc.d13c' %in% names(soc.data))) soc.data[, soc.d13c := init.soc.d13c]
  soc.data[, dpm.d13c := soc.d13c]
  soc.data[, rpm.d13c := soc.d13c]
  soc.data[, doc.d13c := soc.d13c]
  soc.data[, bio.d13c := soc.d13c]
  soc.data[, sta.d13c := soc.d13c]
  soc.data[, co2.d13c := 0.0]
  return(soc.data)
}

#' Initialise Somic data
#'
#' Initialises model inputs
#'
#' @param soc.data A data frame.
#'
#' @return The data.frame but with more stuff filled in.
#' \item{soc.data}{The data.frame but with more stuff filled in.}
#' @author Dominic Woolf
#' @note This is a very simple function.
#' @rdname initialise.daily.data
#' @export
initialise.daily.data <- function(soc.data) {
  # function to intialise values in soc dataframe
  # convert monthly to daily data
  start.date = as.Date('1800-01-01')
  soc.data[, date := start.date %m+% months(time-1)]
  daily.soc.data = soc.data[, .(date = seq(date[1], date[.N], by='1 day')), by = exp]
  daily.soc.data = merge(soc.data, daily.soc.data, by=c('exp', 'date'), all.y = T)
  setorder(daily.soc.data, exp, date)
  
  #interpolate missing daily values
  # daily.soc.data[, month := elapsed_months(date, start.date) + 1]
  daily.soc.data[, time := 1:.N, by = exp]
  daily.soc.data[, date := NULL]
  daily.soc.data[, temp := na.spline(temp), by = exp]
  daily.soc.data[, precip := na.spline(precip)*12/365, by = exp]
  daily.soc.data[, pet := na.spline(pet)*12/365, by = exp]
  daily.soc.data[, cover := na.locf(cover), by = exp]
  daily.soc.data[, clay := na.locf(clay), by = exp]
  daily.soc.data[, depth := na.locf(depth), by = exp]
  daily.soc.data[, soc := na.locf(soc), by = exp]
  daily.soc.data[, soc.d13c := na.locf(soc.d13c), by = exp]
  daily.soc.data[, max_tsmd := na.locf(max_tsmd), by = exp]
  daily.soc.data[, max_tsmd := na.locf(max_tsmd), by = exp]
  daily.soc.data[is.na(added.dpm), added.dpm := 0.0]
  daily.soc.data[is.na(added.rpm), added.rpm := 0.0]
  daily.soc.data[is.na(added.bio), added.bio := 0.0]
  daily.soc.data[is.na(added.hum), added.hum := 0.0]
  daily.soc.data[, added.d13c := na.locf(added.d13c), by = exp]
  daily.soc.data = initialise.somic.data(daily.soc.data)
  daily.soc.data
}

#' Initialise RothC data
#'
#' Initialises model inputs
#'
#' @param soc.data A data frame.
#'
#' @return A data.table but with more stuff filled in.
#' \item{soc.data}{The data.frame but with more stuff filled in.}
#' @author Dominic Woolf
#' @note This is a very simple function to intialise the data structure.
#' @rdname initialise.rothc.data
#' @export
initialise.rothc.data <- function(soc.data) {
  iom = soc.data[1, 0.049 * soc^1.139]
  soc.data[, soc := soc - iom]
  soc.data[, dpm := soc * 0.001]
  soc.data[, rpm := soc * 0.081]
  soc.data[, bio := soc * 0.0331]
  soc.data[, sta := soc * (1-0.001-0.081-0.0331)]
  soc.data[, atsmd := 0.0]
  soc.data[, a := fT.RothC(temp)]
  soc.data[, c := ifelse (cover==0, 1, 0.6)]
  soc.data
}

#' Weighted mean of raster bricks
#'
#' Finds the weighted mean of bricks.
#'
#' @param b.list A list of bricks.
#' @param w A vector of weights- same length as \code{b.list}.
#'
#' @return A ratser brick containing the weighted means.
#' \item{wmean.bricks}{The sum of the squared values.}
#' \item{b.list}{A list of raster bricks.}
#' \item{w}{weights (same length as b.list).}
#' @author Dominic Woolf
#' @note This is a very simple function.
#' @rdname wmean.bricks
#' @export
wmean.bricks = function(b.list, w) {
  b.wmean = b.list[[1]] * w[1]
  for (i in seq_along(w)[-1]){
    b.wmean = b.wmean + b.list[[i]] * w[i]
  }
  b.wmean / sum(w)
}


#' Inverse of cumsum
#'
#' Finds the inverse of cumsum.
#'
#' @param x A vector.
#'
#' @return A vector.
#' \item{cdiff}{The inverse of cumsum.}
#' \item{x}{A vector.}
#' @author Dominic Woolf
#' @rdname cdiff
#' @export
cdiff = function(x) diff(c(0,x))
