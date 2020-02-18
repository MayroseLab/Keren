R = 10
one_over_R = 1 / R
num_samples = 8
# we would like to get 4 samples between one_over_R and 1 (slow down) and 4 smaples between 1 and R (speed up)
# sampling in equal intervals from [one_over_R, R] will give too many speed up samples...
# convert to log, sample uniform, take exponent:
log_R = log(R)
log_one_over_R = log(one_over_R)
step_size = (log_R - log_one_over_R) / (num_samples - 1)
log_samples = seq(log_one_over_R,log_R,by=step_size)
samples = exp(log_samples)
print(samples)
# same factor between two sampled values:
for (i in 2:length(samples))
{
  factor_i = samples[i] / samples[i-1]
  print(factor_i)
}