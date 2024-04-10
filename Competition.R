m = matrix(sample(1:100), nrow = 100, ncol =100)
submatrix_size = 5
highest_mean <- -Inf
highest_mean_submatrix <- NULL

for (r in 1:(nrow(m) - submatrix_size + 1)){
  for (c in 1:(ncol(m) - submatrix_size +1)){
    submatrix = m[r:(r + submatrix_size -1), c:(c + submatrix_size -1)]
  submatrix_mean = mean(submatrix)
  if (submatrix_mean > highest_mean){
    highest_mean = submatrix_mean
    highest_mean_submatrix = submatrix
    }
  }
}
print(highest_mean)
print(highest_mean_submatrix)