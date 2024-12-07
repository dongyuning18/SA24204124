#' @title A illustration dataset
#' @name data
#' @description Generate a dataset, SP500 returns from https://finance.yahoo.com/, GDP growth, Inflation, Unemployment rate of China from https://data.worldbank.org.cn.
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' plot(Year, SP_Return)
#' plot(Year, GDP)
#' plot(Year, Inflation)
#' plot(Year, Unemployment)
#' }
NULL

#' @title Stock risk of data
#' @name dataprediction
#' @description Using yearly sp500 returns to predict the stock risk value in China, by using data.
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' r1 = predict_stock_market_risk(SP_Return,GDP,Inflation,Unemployment,window_size=10)$risk_predictions
#' plot(seq(2024,2038),r1,ylab = 'Stock Risk Value',type = 'l')
#' }
#' @import knitr
#' @import rmarkdown
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm sd coef lm
#' @import TSA
#' @import DAAG
#' @import boot
#' @import MASS
#' @import bootstrap
#' @import cramer
#' @import coda
#' @import nloptr
#' @import grid
#' @import jpeg
#' @import microbenchmark
#' @useDynLib SA24204124
NULL


#' @title Forecast stock market risk based on macroeconomic data
#' @description Combine gdp, unemployment rate and inflation rate to predict future stock market risk
#' @param stock_returns stock returns
#' @param x1 The first factor, such as GDP growth
#' @param x2 The second factor, such as inflation rate
#' @param x3 The third factor, such as umemployment rate
#' @param window_size The size of the scroll window, 50 by default
#' @return prediction of risk value of size n-window_size+1.
#' @examples
#' \dontrun{
#' stock_returns <- rnorm(200)
#' gdp_growth <- rnorm(200)
#' inflation_rate <- rnorm(200)
#' unemployment_rate <- rnorm(200) 
#' predict_stock_market_risk(stock_returns, gdp_growth, inflation_rate, unemployment_rate)
#' }
#' @export
predict_stock_market_risk <- function(stock_returns, x1, x2, x3, window_size = 80) {
  n <- length(stock_returns)
  
  #检查输入的参数长度是否一致
  if (n != length(x1) || n != length(x2) || n != length(x3)) {
    stop("All input data vectors must have the same length.")
  }
  
  # 检查窗口大小是否合法
  if (window_size <= 0 || window_size > n) {
    stop("Window size must be a positive integer and smaller than or equal to the length of the data.")
  }
  
  #存储风险预测结果
  risk_predictions <- numeric(n - window_size + 1)
  
  for (i in 1:(n - window_size + 1)) {
    # 提取当前窗口的数据
    current_window_stock_returns <- stock_returns[i:(i + window_size - 1)]
    current_window1 <- x1[i:(i + window_size - 1)]
    current_window2 <- x2[i:(i + window_size - 1)]
    current_window3 <- x3[i:(i + window_size - 1)]
    
    # 检查窗口内是否包含NA值
    if (any(is.na(current_window_stock_returns)) || any(is.na(current_window1)) || 
        any(is.na(current_window2)) || any(is.na(current_window3))) {
      risk_predictions[i] <- NA  # 如果窗口内有NA值，设置风险为NA
      next
    }
    
    # 计算当前窗口股市收益率的标准差作为风险
    current_risk <- sd(current_window_stock_returns)
    
    # 将经济因素（x1, x2, x3）合并成一个数据框
    data <- data.frame(current_risk = rep(current_risk, length(current_window1)), x1 = current_window1, x2 = current_window2, x3 = current_window3)
    
    
    # 多元线性回归：股市风险 ~ GDP + 通货膨胀率 + 失业率
    model <- lm(current_risk ~ x1 + x2 + x3, data = data)
    
    # 提取回归模型的系数
    coefficients <- coef(model)
    
    # 计算回归模型的预测值（即股市风险）
    predicted_risk <- sum(coefficients[-1] * c(mean(current_window1), mean(current_window2), mean(current_window3))) + coefficients[1]
    
    # 存储预测的风险值
    risk_predictions[i] <- predicted_risk
  }
  
  # 输出预测结果和统计信息
  result <- list(
    risk_predictions = risk_predictions,
    mean_risk = mean(risk_predictions, na.rm = TRUE),
    sd_risk = sd(risk_predictions, na.rm = TRUE),
    min_risk = min(risk_predictions, na.rm = TRUE),
    max_risk = max(risk_predictions, na.rm = TRUE)
  )
  
  return(result)
}


#' @title An illustration dataset with 2 columns: pressure and temp to be detected
#' @name df
#' @description A simple example dataset, detect the outliers of first column temp and the 2nd column pressure.
#' @examples
#' \dontrun{
#' data(df)
#' attach(df)
#' plot(temp)
#' plot(pressure)
#' }
NULL

#' @title Given range of dataset df.
#' @name ranges
#' @description Two lists, each llist is a binary vector representing range.
#' @examples
#' \dontrun{
#' data(ranges)
#' attach(ranges)
#' print(ranges)
#' }
NULL
