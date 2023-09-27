
repVec = function(vec, nums) {
  n = length(vec)
  m = length(nums)
  if (n != m) {
    stop("Argument lengths don't match")
  }
  resLength = sum(nums)
  res = rep(vec[1], resLength)
  start = 0
  for (i in 1:n) {
    end = start + nums[i]
    res[(start+1):end] = rep(vec[i], nums[i])
    start = end
  }
  return(res)
}

repMatrix = function(matrix, nums) {
  n = dim(matrix)[1]
  m = length(nums)
  if (n != m) {
    stop("Argument lengths don't match")
  }
  resLength = sum(nums)
  res = matrix(rep(0, resLength*dim(matrix)[2]), nrow = resLength)
  start = 0
  for (i in 1:n) {
    end = start + nums[i]
    res[(start+1):end,] = t(matrix(rep(matrix[i,], nums[i]), ncol = nums[i]))
    start = end
  }
  return(res)
}

sumVec = function(vec, nums) {
  n = length(vec)
  if (n != sum(nums)) {
    stop("Argument length does not match window lengths")
  }
  m = length(nums)
  res = rep(0, m)
  start = 0
  for (i in 1:m) {
    end = start+nums[i]
    res[i] = sum(vec[(start+1):end])
    start = end
  }
  return(res)
}