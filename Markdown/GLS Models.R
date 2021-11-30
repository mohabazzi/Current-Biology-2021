gls.model.full <- function(data, corStruct = FALSE, acf.plot = FALSE, method = "ML", p.adjust.type = "none", rando = FALSE) {
  if(acf.plot) par(mfrow = c(3,2), mar = c(2,4,3,1))
  
  # 1. Null Model
  null <- gls(N~1, data = data, method = method)
  if(corStruct) {
    corr <- Initialize(corAR1(form = ~time), data = data)
    null <- gls(N~1, data = data, correlation = corr, method = method)
  }
  null.resid <- residuals(null, type = "normalized")
  null.rhos <- acf(null.resid, plot = acf.plot, lag.max = 12, main = "Null")
  null.LBtest <- Box.test(null.resid, lag=12, type="Ljung-Box", fitdf = 1)
  null.sum <- summary(null)
  
  # 2. Sea-level
  sl <- gls(N~sea.level, data = data, method = method)
  if(corStruct) {
    corr <- Initialize(corAR1(form = ~time), data = data)
    sl <- gls(N~sea.level, data = data, correlation = corr, method = method)
  }
  sl.resid <- residuals(sl, type = "normalized")
  sl.rhos <- acf(sl.resid, plot = acf.plot, lag.max = 12, main = "Sea-Level")
  sl.LBtest <- Box.test(sl.resid, lag=12, type="Ljung-Box", fitdf = 1)
  sl.sum <- summary(sl)
  
  # 3. Temperature
  temp <- gls(N~delta.O18, data = data, method = method)
  if(corStruct) {
    corr <- Initialize(corAR1(form = ~time), data = data)
    temp <- gls(N~delta.O18, data = data, correlation = corr, method = method)
  }
  temp.resid <- residuals(temp, type = "normalized")
  temp.rhos <- acf(temp.resid, plot = acf.plot, lag.max = 12, main = "Temperature")
  temp.LBtest <- Box.test(temp.resid, lag=12, type="Ljung-Box", fitdf = 1)
  temp.sum <- summary(temp)
  
  # 4. Sea Level + Temperature
  sltemp <- gls(N~sea.level + delta.O18, data = data, method = method)
  if(corStruct) {
    corr <- Initialize(corAR1(form = ~time), data = data)
    sltemp <- gls(N~sea.level + delta.O18, data = data, correlation = corr, method = method)
  }
  sltemp.resid <- residuals(sltemp, type = "normalized")
  sltemp.rhos <- acf(sltemp.resid, plot = acf.plot, lag.max = 12, main = "Sea-Level + Temperature")
  sltemp.LBtest <- Box.test(sltemp.resid, lag=12, type="Ljung-Box", fitdf = 1)
  sltemp.sum <- summary(sltemp)
  
  # 5. Sea Level-Temperature Interaction
  sltemp.int <- gls(N~sea.level * delta.O18, data = data, method = method)
  if(corStruct) {
    corr <- Initialize(corAR1(form = ~time), data = data)
    sltemp.int <- gls(N~sea.level * delta.O18, data = data, correlation = corr, method = method)
  }
  sltemp.int.resid <- residuals(sltemp.int, type = "normalized")
  sltemp.int.rhos <- acf(sltemp.int.resid, plot = acf.plot, lag.max = 12, main = "Sea-Level * Temperature")
  sltemp.int.LBtest <- Box.test(sltemp.int.resid, lag=12, type="Ljung-Box", fitdf = 1)
  sltemp.int.sum <- summary(sltemp.int)
  
  # Results
  mods <- c("Null", "Sea-level", "Delta.18O", "Sea-level+Delta.18O","Sea-level*Delta.18O")
  acf.res <- list(null.rhos, sl.rhos, temp.rhos, sltemp.rhos, sltemp.int.rhos)
  names(acf.res) <- mods
  lag1.res <- data.frame(lag1.acf = c(null.rhos$acf[2], sl.rhos$acf[2], temp.rhos$acf[2], sltemp.rhos$acf[2],sltemp.int.rhos$acf[2]),
                         pvalue = p.adjust(c(null.LBtest$p.value, sl.LBtest$p.value, temp.LBtest$p.value,sltemp.LBtest$p.value,sltemp.int.LBtest$p.value), method = p.adjust.type),
                         row.names = mods)
  aics <- c(null.sum$AIC, sl.sum$AIC, temp.sum$AIC,sltemp.sum$AIC,sltemp.int.sum$AIC)
  aicc <- c(MuMIn::AICc(null), MuMIn::AICc(sl), MuMIn::AICc(temp),MuMIn::AICc(sltemp),MuMIn::AICc(sltemp.int))
  ord <- order(aicc)
  mod.list <- list(null, sl, temp, sltemp, sltemp.int)[ord]
  lik.rat <- anova(mod.list[[1]], mod.list[[2]], mod.list[[3]], mod.list[[4]], mod.list[[5]])
  lik.res <- data.frame(lik.rat[,-1], 
                        AIC.W = round(akaike.weights(aics)$weights, 4)[ord], 
                        AICc = aicc[ord], 
                        AICc.W = round(akaike.weights(aicc)$weights, 4)[ord], Model.Factor = factor(mods[ord]))
  res.list <- list(ACF.res = acf.res[ord], LAG1.res = lag1.res,
                   gls.models = list(Null = null, 
                                     'Sea-Level' = sl, 
                                     Temperature = temp,
                                     'SL+Temp' = sltemp,
                                     'SL*Temp' = sltemp.int),
                   gls.coefficients = list(Null = coef(null.sum), 
                                           'Sea-Level' = coef(sl.sum), 
                                           Temperature = coef(temp.sum),
                                           'SL+Temp' = coef(sltemp.sum),
                                           'SL*Temp' = coef(sltemp.int.sum)),
                   Overall.Results = lik.res)
}