#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stats)
  library(pROC)
})

# 1. Setup Mock Data
set.seed(42)
n <- 10000
df <- data.frame(
  x = rnorm(n),
  y = rnorm(n) * 2 + 5,
  cat1 = sample(c('A', 'B', 'C'), n, replace = TRUE),
  cat2 = sample(c('X', 'Y'), n, replace = TRUE),
  binary = rbinom(n, 1, 0.5)
)
df_missing <- df
df_missing$x[10:20] <- NA

# merge needs a key column: Stats::LikeR has no row.names to join on, so all
# three scripts join two copies of the frame on an explicit id.  Joining by
# row.names, as this script used to, is about half the work of a real key join.
df_id <- df
df_id$id <- seq_len(n)

y_true <- rbinom(n, 1, 0.5)
y_scores <- runif(n)

# 2. Define Function Mappings (Stats::LikeR -> R equivalent)
benchmarks <- list(
  # Basic Stats
  list(liker="mean", r_sub="base::mean", func=function() mean(df$x)),
  list(liker="median", r_sub="stats::median", func=function() median(df$x)),
  list(liker="sd", r_sub="stats::sd", func=function() sd(df$x)),
  list(liker="var", r_sub="stats::var", func=function() var(df$x)),
  list(liker="min", r_sub="base::min", func=function() min(df$x)),
  list(liker="max", r_sub="base::max", func=function() max(df$x)),
  list(liker="quantile", r_sub="stats::quantile", func=function() quantile(df$x, probs=c(0.25, 0.5, 0.75))),
  # one coefficient, not a 2x2 matrix: Stats::LikeR's cor/cov take two vectors
  list(liker="cor", r_sub="stats::cor", func=function() cor(df$x, df$y)),
  list(liker="cov", r_sub="stats::cov", func=function() cov(df$x, df$y)),
  
  # Distributions & Tests
  list(liker="rnorm", r_sub="stats::rnorm", func=function() rnorm(1000)),
  list(liker="runif", r_sub="stats::runif", func=function() runif(1000)),
  list(liker="t_test", r_sub="stats::t.test", func=function() t.test(x ~ cat2, data = df)),
  list(liker="wilcox_test", r_sub="stats::wilcox.test", func=function() wilcox.test(x ~ cat2, data = df)),
  list(liker="chisq_test", r_sub="stats::chisq.test", func=function() chisq.test(table(df$cat1, df$cat2))),
  list(liker="shapiro_test", r_sub="stats::shapiro.test", func=function() shapiro.test(df$x[1:5000])),
  list(liker="binom_test", r_sub="stats::binom.test", func=function() binom.test(500, 1000, 0.5)),
  
  # Data Manipulation
  list(liker="filter", r_sub="dplyr::filter", func=function() filter(df, x > 0)),
  list(liker="select_cols", r_sub="dplyr::select", func=function() select(df, x, cat1)),
  list(liker="drop_cols", r_sub="dplyr::select", func=function() select(df, -y)),
  list(liker="rename_cols", r_sub="dplyr::rename", func=function() rename(df, Category_1 = cat1)),
  list(liker="dropna", r_sub="stats::na.omit", func=function() na.omit(df_missing)),
  # every column, as df.fillna(0) and fillna($df, value => 0) both do
  list(liker="fillna", r_sub="tidyr::replace_na", func=function() mutate(df_missing, across(everything(), ~ replace_na(.x, 0)))),
  list(liker="drop_duplicates", r_sub="dplyr::distinct", func=function() distinct(df)),
  list(liker="group_by", r_sub="dplyr::group_by", func=function() df %>% group_by(cat1) %>% summarise(mean_x = mean(x))),
  list(liker="concat", r_sub="base::rbind", func=function() rbind(df, df)),
  list(liker="merge", r_sub="base::merge", func=function() merge(df_id, df_id, by="id")),
  list(liker="value_counts", r_sub="base::table", func=function() table(df$cat1)),
  list(liker="pivot_table", r_sub="tidyr::pivot_wider", func=function() df %>% group_by(cat1, cat2) %>% summarise(val=mean(x), .groups='drop') %>% pivot_wider(names_from=cat2, values_from=val)),
  
  # Modeling & Metrics
  list(liker="lm", r_sub="stats::lm", func=function() lm(y ~ x + cat1, data=df)),
  list(liker="glm", r_sub="stats::glm", func=function() glm(binary ~ x + y, family=binomial, data=df)),
  list(liker="auc", r_sub="pROC::auc", func=function() auc(y_true, y_scores, quiet=TRUE)),
  list(liker="prcomp", r_sub="stats::prcomp", func=function() prcomp(df[, c("x", "y")]))
)

# 3. Execution Engine
runs <- 7
results_list <- list()

cat("Running R benchmarks...\n")
for (b in benchmarks) {
  for (i in 1:runs) {
    # Force garbage collection, then record the baseline heap so the reported
    # figure is the allocation caused by this operation rather than the
    # session-wide high-water mark (which would include df itself).
    # Ncells are 56 bytes and Vcells 8 bytes on 64-bit builds.
    gc(reset = TRUE, full = TRUE)
    base_mem <- sum(gc(full = TRUE)[, "used"] * c(56, 8))

    start_time <- Sys.time()

    tryCatch({
      res <- b$func()
    }, error = function(e) {
      cat("Error in", b$liker, ":", conditionMessage(e), "\n")
    })

    end_time <- Sys.time()

    # Peak heap since the reset above, minus baseline, in bytes. Coarser than
    # Python's tracemalloc (R's GC reports at page granularity), but the same
    # order of magnitude rather than ~1000x off.
    peak_mem <- sum(gc(full = TRUE)[, 6] * 1024^2) - base_mem

    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    results_list[[length(results_list) + 1]] <- data.frame(
      `Stats::LikeR function` = b$liker,
      `Language equivalent` = b$r_sub,
      time = elapsed,
      RAM = peak_mem,
      check.names = FALSE
    )
  }
}

# 4. Output
final_results <- bind_rows(results_list)
output_file <- "r_benchmarks.tsv"
write.table(final_results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("Done. R results written to %s\n", output_file))