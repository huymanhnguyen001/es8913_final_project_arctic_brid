# Function for replacing zero with NA values for plotting scatter plot
zero_it_out <- function (data, col, options) {
  if (options == "by_NA") {
    for (col in col) {
      data[,col][data[,col] == 0] <- NA
    }
    return(data)
  }
}
