# Regenerates the frozen table in t/cor_test.pearson.ci.R.t.
#   Rscript t/cor_test.pearson.ci.R.R
# The test never runs this; it reads the literals that were pasted from it.
#
# Emits one row per (dataset, alternative, conf.level): the Fisher-z interval
# R reports, plus the statistic and p-value alongside it.  A one-sided
# alternative gets a one-sided interval in R -- tanh(-Inf) and tanh(Inf), i.e.
# exactly -1 and 1 for the open end -- which is the whole point of the table.
options(digits = 17, scipen = 500)

cases <- list()
add <- function(name, x, y) cases[[length(cases)+1]] <<- list(name=name, x=x, y=y)
add("strong10", 1:10, c(2,1,4,3,7,5,9,6,10,8))
add("weak12",   1:12, c(5,2,9,1,12,3,11,4,8,6,10,7))
add("neg8",     1:8,  c(8,6,7,5,3,4,1,2))
add("n4",       1:4,  c(2,1,4,3))
add("n5flat",   1:5,  c(1,3,2,5,4))
set.seed(4)
add("noisy40",  1:40, round(rnorm(40), 6))

for (cs in cases) {
  for (alt in c("two.sided", "greater", "less")) {
    for (cl in c(0.9, 0.95, 0.99)) {
      r <- cor.test(cs$x, cs$y, method="pearson", alternative=alt, conf.level=cl)
      cat(sprintf("%s\t%s\t%s\t%s\t%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
                  cs$name, paste(cs$x, collapse=","), paste(cs$y, collapse=","),
                  alt, cl, r$statistic, r$p.value, r$estimate,
                  r$conf.int[1], r$conf.int[2]))
    }
  }
}
