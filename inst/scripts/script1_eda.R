### exploratory data analysis 1

alt <- coords$altitude[order(coords$site)]

# mean
png("img/mean_Tmax_by_site.png", width = 400, height = 300)
y <- tapply(df0$Tmax, df0$site, mean)
plot(x = alt, y,
     xlab = "Elevation", ylab = "Mean Tmax by site")
abline(lm(y ~ alt))
dev.off()

# sd
png("img/sd_Tmax_by_site.png", width = 400, height = 300)
y <- tapply(df0$Tmax, df0$site, sd)
plot(x = alt, y,
     xlab = "Elevation", ylab = "sd Tmax by site")
abline(lm(y ~ alt))
dev.off()

# seasonality
png("img/season_Tmax_by_site_yday.png", width = 400, height = 300)
y <- tapply(df0$Tmax, list(df0$yday, df0$site), mean)
matplot(
  x = as.numeric(rownames(y)),
  y = y,                      
  type = "l",                 
  lty = 1,                    
  col = 1:ncol(y),            
  xlab = "Day within year",
  ylab = "Mean Tmax by site and yday"
)
dev.off()

# long-term trend
png("img/mean_Tmax_by_site_year.png", width = 400, height = 300)
y <- tapply(df0$Tmax, list(df0$year, df0$site), mean)
matplot(
  x = as.numeric(rownames(y)),
  y = y,                      
  type = "l",                 
  lty = 1,                    
  col = 1:ncol(y),            
  xlab = "Year",
  ylab = "Mean Tmax by site and year"
)
dev.off()

# autoregression
png("img/ar_Tmax_by_site.png", width = 400, height = 300)
y <- tapply(df0$Tmax, df0$site, function(x) acf(x)$acf[,,1][2])
plot(x = alt, y,
     xlab = "Elevation", ylab = "Autoregression Tmax by site")
dev.off()



### exploratory data analysis 2

tau <- c(0.05, 1:9/10, 0.95)

# quantile
quantiles <- tapply(df0$Tmax,
                    df0$site,
                    quantile, probs = tau)

xxx <- t(sapply(quantiles, function(x) x))

png("img/q_Tmax.png", width = 400, height = 300)
boxplot(xxx, xlab = "Quantile level (tau)", ylab = "Quantile Tmax")
dev.off()

# quantile long-term trend
quantiles_06_15 <- tapply(df0$Tmax[df0$year %in% 2006:2015],
                          df0$site[df0$year %in% 2006:2015],
                          quantile, probs = tau)

quantiles_76_85 <- tapply(df0$Tmax[df0$year %in% 1976:1985],
                          df0$site[df0$year %in% 1976:1985],
                          quantile, probs = tau)

xxx <- t(mapply(function(x, y) x - y, quantiles_06_15, quantiles_76_85))

png("img/q_Tmax_by_site_year.png", width = 400, height = 300)
boxplot(xxx, xlab = "Quantile level (tau)", ylab = "Difference Tmax")
dev.off()
