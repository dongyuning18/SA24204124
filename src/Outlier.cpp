#include <Rcpp.h>
using namespace Rcpp;

//' @title rangeAnomalyDetection
//' @description Outlier detection by range
//' @param df data to be detected if there are outliers, type of data frame
//' @param value_ranges data out of this normal range will be extracted, IQR(interquantile range) by default.
//' @return a list, including the index of anomaly, and the range.
//' @examples
//' \dontrun{
//' df2 <- data.frame(
//' d1 = runif(10, -30, 40),
//' d2 = runif(10, 900, 1100)
//' )
//' rangeAnomalyDetection(df2, value_ranges)
//' }
//' @export
// [[Rcpp::export]]
List rangeAnomalyDetection(DataFrame df, Nullable<List> value_ranges = R_NilValue) {
  
  // 创建一个空的列表来存储异常值
  List anomalies;
  
  // 如果没有提供范围，自动基于 IQR 计算范围
  if (value_ranges.isNull()) {
    List ranges;
    
    // 获取 quantile 函数
    Function quantile("quantile");
    
    for (int i = 0; i < df.size(); i++) {
      CharacterVector colnames = df.names();
      String colname = colnames[i];
      NumericVector column = df[colname];
      
      if (Rf_isNumeric(column)) {
        // 调用 R 中的 quantile 函数计算 Q1 和 Q3
        NumericVector Q1 = quantile(column, 0.25);
        NumericVector Q3 = quantile(column, 0.75);
        double IQR = Q3[0] - Q1[0];
        
        double lower_bound = Q1[0] - 1.5 * IQR;
        double upper_bound = Q3[0] + 1.5 * IQR;
        
        NumericVector range = NumericVector::create(lower_bound, upper_bound);
        ranges.push_back(range, colname);
      }
    }
    
    // 将 ranges 转换为 Nullable 类型
    value_ranges = wrap(ranges);
  }
  
  // 获取列名
  CharacterVector colnames = df.names();
  
  // 标记超出范围的异常值
  for (int i = 0; i < df.size(); i++) {
    String colname = colnames[i];
    NumericVector column = df[colname];
    
    if (Rf_isNumeric(column)) {
      // 使用 as<List>() 解包 value_ranges，获取列表
      List unwrapped_value_ranges = as<List>(value_ranges);
      
      // 获取指定列的范围
      NumericVector range = unwrapped_value_ranges[colname];
      double lower_bound = range[0];
      double upper_bound = range[1];
      
      // 手动实现 which()，找出异常值的索引
      IntegerVector outlier_indices;
      for (int j = 0; j < column.size(); j++) {
        if (column[j] < lower_bound || column[j] > upper_bound) {
          outlier_indices.push_back(j + 1);  // 索引从 1 开始
        }
      }
      
      anomalies.push_back(outlier_indices, colname);
    }
  }
  
  // 返回包含异常值的记录和计算的范围
  return List::create(Named("anomalies") = anomalies,
                      Named("value_ranges") = value_ranges);
}
