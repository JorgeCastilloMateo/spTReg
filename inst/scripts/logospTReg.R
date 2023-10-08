# SCRIPT TO REPRODUCE spTReg LOGO (R VERSION 4.3.1)
set.seed(23)
df <- data.frame(x = stats::rnorm(1000, 0, 0.5))
#df$y <- 10 + 4 * df$x^2 + stats::rnorm(1000, sd = 15)
df$y <- df$x * sin(1 / df$x) + stats::rnorm(1000, sd = 0.1)

for (tau in 1:3) {
  df[, tau + 2] <- colMeans(spTReg::iidm(y ~ I(x * sin(1 / x)), data = df, 
    method = "quantile", quantile = c(0.01, 0.50, 0.99)[tau],
    verbose = FALSE)$fitted.values)
}

p <- ggplot2::ggplot(data = df) + ggplot2::xlim(max(abs(df$x)) * c(-1.01, 1.01)) +
  ggplot2::geom_point(ggplot2::aes(x = x, y = y), size = 0.2, alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(x = x, y = V3), size = 0.2, col = "red") +
  ggplot2::geom_line(ggplot2::aes(x = x, y = V4), size = 0.2, col = "red") +
  ggplot2::geom_line(ggplot2::aes(x = x, y = V5), size = 0.2, col = "red")

p <- p + 
  ggplot2::labs(
    title    = ggplot2::element_blank(), 
    subtitle = ggplot2::element_blank()) +
  ggplot2::theme(
    axis.line        = ggplot2::element_blank(), 
    axis.text.x      = ggplot2::element_blank(),
    axis.text.y      = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.title.x     = ggplot2::element_blank(),
    axis.title.y     = ggplot2::element_blank(),
    legend.position  = "none",
    panel.background = ggplot2::element_blank(),
    panel.border     = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.background  = ggplot2::element_blank())

if (!require("hexSticker")) install.packages("hexSticker")
hexSticker::sticker(p, package = "spTReg", 
  s_x = 1, s_y = 0.8, s_width = 1.7, s_height = 1.35,
  p_x = 1, p_y = 1.4, p_size = 23, h_size = 1.5,
  h_fill = "grey", 
  h_color = "black",
  filename = paste0(getwd(), "/inst/img/logospTReg.png"))
