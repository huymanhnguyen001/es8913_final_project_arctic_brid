zero_it_out <- function (data, col, options) {
  if (options == "by_NA") {
    for (col in col) {
      data[,col][data[,col] == 0] <- NA
    }
    return(data)
  }
  if (options == "LOD") {
    data <- as.data.frame(data)
    for (col in col) {
      data[,col][data[,col] == 0] <- runif(sum(data[,col] == 0),
                                                     min = 0.0145000,
                                                     max = 0.0155000)
    }
    return(data)
  }
  # if (options == "no_zero") {
  #   list <- list()
  #   for (col in col) {
  #     list[col] <- filter(data[,col] > 0)
  #   }
  #   
  # }
}
