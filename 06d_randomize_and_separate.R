library(kgp)
data(kgp)

kgp3


# Convert to data.table if not already
dt <- as.data.table(kgp3)

# Take one sample per combination of sex, pop, and reg
representative_sample <- dt[, .SD[sample(.N, 1)], by = .(sex, pop, reg)]

# View the result
print(representative_sample$id)
