#!/usr/bin/env Rscript
# The R side of plot.scaling.pl: the same functions, the same size ladder.
#
# Where benchmark.R asks "how long does each function take on one 10,000-row
# frame", this asks "what shape is that curve".  Each call is timed at a ladder
# of sizes half a decade apart, seven times at each size, and written out for
# plot.scaling.pl to draw against the same measurements from Perl and Python.
# On a log-log axis the slope of the line is the exponent.
#
#     perl plot.scaling.pl --data                 # once, writes the fixtures
#     perl -Iblib/arch -Iblib/lib \
#          plot.scaling.pl --measure              # -> perl_scaling.tsv
#     python3 scale.py                            # -> python_scaling.tsv
#     Rscript scale.R                             # -> r_scaling.tsv
#     perl plot.scaling.pl --plot                 # -> scaling.*.svg
#
# Environment:
#
#     SCALE_DIR    where --data put the fixtures (/tmp/likeR.scaling)
#     SCALE_RUNS   runs per (function, size); default 7
#     SCALE_CAP    seconds; once one run of a function takes longer than this,
#                  that function is not tried at any larger size.  Default 4.
#     SCALE_MAX_N  hard ceiling on the row count, for a quick partial run
#     SCALE_TARGET seconds a single measurement should span; a call faster than
#                  this is repeated until it does.  Default 0.002.
#     SCALE_CPU    the CPU this run pins itself to; default 0, "none" to let
#                  the scheduler place it.
#
# On a hybrid CPU -- Intel's P-core/E-core parts, and the big.LITTLE ARM
# designs -- a process landing on an E-core reads 1.5 to 2 times slower than
# the identical loop on a P-core, which is a wider band than most of the
# differences these plots are drawn to show.  So this pins itself to one CPU
# rather than asking to be started under "taskset"; see pin_to_one_cpu() below.
# plot.scaling.pl --measure and scale.py do the same, and all three default to
# CPU 0, so the three languages land on the same core without anyone having to
# remember.  SCALE_CPU moves it, and has to be given the same value in all
# three.
#
# Four panels have no R line, and are left empty rather than filled with
# something else:
#
#   skew, kurtosis   base R has neither.  e1071 and moments both provide them,
#                    with three different type= conventions between them, so
#                    there is no single call that is "R's answer".
#   aoh2hoa          there is no row-record frame in R to convert from.  The
#                    nearest thing, do.call(rbind, lapply(rows, as.data.frame)),
#                    is quadratic and would be a comment on that idiom rather
#                    than on R.
#
# Two panels do share a call with another panel, deliberately, and it is noted
# at each: R is columnar, so read.csv is the honest counterpart of both
# read_table's default row-record output and its output.type => 'hoa', and the
# difference between those two panels is the cost of the Perl shape alone.
suppressPackageStartupMessages({
    library(dplyr)
    library(pROC)
})

DIR <- Sys.getenv("SCALE_DIR", "/tmp/likeR.scaling")
RUNS <- as.integer(Sys.getenv("SCALE_RUNS", "7"))
CAP <- as.numeric(Sys.getenv("SCALE_CAP", "4"))
MAX_N <- as.integer(Sys.getenv("SCALE_MAX_N", "0"))
TARGET <- as.numeric(Sys.getenv("SCALE_TARGET", "0.002"))
MAX_REPS <- 10000L

# Size ladders.  Half-decade steps; plot.scaling.pl and scale.py carry the same
# three lists, and the fixture row counts are IO_N.  Change one, change all three.
VEC_N <- c(1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6)
IO_N <- c(1e3, 3e3, 1e4, 3e4, 1e5, 3e5)
FRAME_N <- c(1e3, 3e3, 1e4, 3e4, 1e5, 3e5)

if (!is.na(MAX_N) && MAX_N > 0) {
    VEC_N <- VEC_N[VEC_N <= MAX_N]
    IO_N <- IO_N[IO_N <= MAX_N]
    FRAME_N <- FRAME_N[FRAME_N <= MAX_N]
}

# ---------------------------------------------------------------------------
# Input builders -- one per figure, called once per size
# ---------------------------------------------------------------------------
build_vector <- function(n) {
    set.seed(42)
    list(x = rnorm(n), y = rnorm(n, 5, 2),
         label = rbinom(n, 1, 0.5), n = n)
}

build_io <- function(n) {
    f <- list(
        num_csv = file.path(DIR, sprintf("num.%d.csv", n)),
        mix_csv = file.path(DIR, sprintf("mix.%d.csv", n)),
        mix_tsv = file.path(DIR, sprintf("mix.%d.tsv", n)),
        out = file.path(DIR, sprintf("out.r.%d.tmp", Sys.getpid()))
    )
    for (k in c("num_csv", "mix_csv", "mix_tsv")) {
        if (!file.exists(f[[k]])) {
            stop(sprintf('missing fixture "%s"; run "perl plot.scaling.pl --data" first',
                         f[[k]]))
        }
    }
    # The frame handed to write.csv is read in, not synthesized, so all three
    # languages write out the same table.
    f$df <- read.csv(f$mix_csv, stringsAsFactors = FALSE)
    f
}

build_frame <- function(n) {
    set.seed(42)
    df <- data.frame(
        x = rnorm(n),
        y = rnorm(n, 5, 2),
        cat1 = sample(c("A", "B", "C"), n, replace = TRUE),
        cat2 = sample(c("X", "Y"), n, replace = TRUE),
        binary = rbinom(n, 1, 0.5),
        stringsAsFactors = FALSE
    )
    # merge needs a key column: a Stats::LikeR frame has no row.names to join
    # on, so all three scripts join on an explicit id, as benchmark.R does.
    df_id <- df
    df_id$id <- seq_len(n)
    list(df = df, df_id = df_id,
         m = as.matrix(df[, c("x", "y", "binary")]), n = n)
}

BUILD <- list(vector = build_vector, io = build_io, frame = build_frame)

# ---------------------------------------------------------------------------
# The benchmarks
# ---------------------------------------------------------------------------
# name is the join key with plot.scaling.pl and scale.py; call is only recorded for
# the reader.  invisible() everywhere, so nothing prints and no print method is
# accidentally inside the timed region.
b <- function(figure, name, call, body) {
    list(figure = figure, name = name, call = call, body = body)
}

BENCHMARKS <- list(
    # --- reductions over one numeric vector --------------------------------
    b("vector", "sum", "sum(x)", function(d) invisible(sum(d$x))),
    b("vector", "min", "min(x)", function(d) invisible(min(d$x))),
    b("vector", "max", "max(x)", function(d) invisible(max(d$x))),
    b("vector", "mean", "mean(x)", function(d) invisible(mean(d$x))),
    b("vector", "median", "median(x)", function(d) invisible(median(d$x))),
    b("vector", "sd", "sd(x)", function(d) invisible(sd(d$x))),
    b("vector", "var", "var(x)", function(d) invisible(var(d$x))),
    b("vector", "quantile", "quantile(x, c(.25,.5,.75))",
      function(d) invisible(quantile(d$x, c(0.25, 0.5, 0.75)))),
    # one number, not the 2x2 matrix cor(cbind(x, y)) would build
    b("vector", "cor", "cor(x, y)", function(d) invisible(cor(d$x, d$y))),
    b("vector", "cov", "cov(x, y)", function(d) invisible(cov(d$x, d$y))),

    # --- transforms that return something the size of their input ----------
    b("transform", "rank", "rank(x)", function(d) invisible(rank(d$x))),
    b("transform", "uniq", "unique(x)", function(d) invisible(unique(d$x))),
    b("transform", "scale", "scale(x)", function(d) invisible(scale(d$x))),
    b("transform", "sample", "sample(x, n %/% 10)",
      function(d) invisible(sample(d$x, d$n %/% 10 + 1))),
    # This one comes out as a flat line, and it is not a measurement error:
    # seq(1, n) with an integer stride returns an ALTREP compact sequence, a
    # pair of numbers that materializes only if something asks it to.  R is
    # being timed on not building the vector at all, where seq() and
    # np.arange() both build it.  Left in for exactly that reason.
    b("transform", "seq", "seq(1, n)", function(d) invisible(seq(1, d$n))),
    b("transform", "auc", "pROC::auc(roc(label, y))",
      function(d) invisible(auc(roc(d$label, d$y, quiet = TRUE)))),

    # --- read_table and write_table, over four inputs each -----------------
    # quote = FALSE throughout: write.csv quotes every character column by
    # default, where write_table and to_csv quote only a field that needs it.
    # None of the fixture's fields do, so quoting would hand R a job the other
    # two were not given.
    b("io", "read_table (csv, numeric)", "read.csv(num.csv)",
      function(d) invisible(read.csv(d$num_csv))),
    b("io", "read_table (csv, mixed)", "read.csv(mix.csv)",
      function(d) invisible(read.csv(d$mix_csv, stringsAsFactors = FALSE))),
    b("io", "read_table (tsv, mixed)", "read.delim(mix.tsv)",
      function(d) invisible(read.delim(d$mix_tsv, stringsAsFactors = FALSE))),
    # the same call as the row-record panel above: R has one table type, so
    # this panel isolates what the Perl shape costs, not what R does differently
    b("io", "read_table (csv, hoa)", "read.csv(mix.csv)",
      function(d) invisible(read.csv(d$mix_csv, stringsAsFactors = FALSE))),
    b("io", "write_table (csv, hoa)", "write.csv(df, f, row.names = FALSE)",
      function(d) write.csv(d$df, d$out, row.names = FALSE, quote = FALSE)),
    b("io", "write_table (tsv, hoa)", 'write.table(df, f, sep = "\\t")',
      function(d) write.table(d$df, d$out, sep = "\t", row.names = FALSE,
                              quote = FALSE)),
    # R has no second table type to write from, so this is the same call again
    b("io", "write_table (csv, aoa)", "write.csv(df, f, row.names = FALSE)",
      function(d) write.csv(d$df, d$out, row.names = FALSE, quote = FALSE)),
    b("io", "write_table (csv, row.names)", "write.csv(df, f, row.names = TRUE)",
      function(d) write.csv(d$df, d$out, row.names = TRUE, quote = FALSE)),

    # --- whole-frame operations --------------------------------------------
    b("frame", "filter", "df[df$x > 0, ]",
      function(d) invisible(d$df[d$df$x > 0, ])),
    b("frame", "select_cols", 'df[, c("x", "cat1")]',
      function(d) invisible(d$df[, c("x", "cat1")])),
    # one column averaged per group, not every numeric one
    b("frame", "group_by + mean", "group_by(cat1) %>% summarise(mean(x))",
      function(d) invisible(d$df %>% group_by(cat1) %>%
                            summarise(m = mean(x), .groups = "drop"))),
    b("frame", "merge", 'merge(df_id, df_id, by = "id")',
      function(d) invisible(merge(d$df_id, d$df_id, by = "id"))),
    b("frame", "value_counts", "table(df$cat1)",
      function(d) invisible(table(d$df$cat1))),
    b("frame", "drop_duplicates", "unique(df)",
      function(d) invisible(unique(d$df))),
    b("frame", "transpose", "t(m)", function(d) invisible(t(d$m)))
)

LADDER <- list(vector = VEC_N, transform = VEC_N, io = IO_N, frame = FRAME_N)
BUILDER_FOR <- list(vector = "vector", transform = "vector", io = "io",
                    frame = "frame")

# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------
# Timed in-process, as scale.py is.  A call faster than TARGET is repeated
# until the pair of clock readings spans TARGET and the total divided by the
# count: Sys.time() resolves microseconds, and min() on a thousand doubles does
# not take one.  The repeat count comes from a second untimed call, so all
# three scripts average over the same span of work.
measure <- function(body, data) {
    t0 <- as.numeric(Sys.time())
    body(data)
    one <- as.numeric(Sys.time()) - t0
    reps <- if (one > 0) min(MAX_REPS, as.integer(TARGET / one) + 1L) else MAX_REPS
    t0 <- as.numeric(Sys.time())
    for (i in seq_len(reps)) body(data)
    list(seconds = (as.numeric(Sys.time()) - t0) / reps, reps = reps)
}

# The CPU this run is pinned to, fixed before the first measurement.
#
# Saying "run this under taskset" in a comment is not the same as running it
# under taskset.  An unpinned rerun of an identical build read Stats::LikeR's
# min() at 2.75 ms against 1.43 pinned, and mean() at 2.61 against 1.35 -- a
# factor of 1.9 on a 13th-generation Core i7, landing on whichever functions
# happened to draw an E-core rather than on all of them evenly.  A reading that
# moves by 1.9x depending on where the scheduler put the process is not a
# measurement, and a cross-language plot drawn from three of them unpinned is
# comparing schedulers.
#
# R has no affinity API -- sched_setaffinity() is not wrapped anywhere in base
# or in parallel -- so this goes through taskset(1) against its own pid, which
# is what plot.scaling.pl does for the same reason.  scale.py needs neither,
# because os.sched_setaffinity is in Python's standard library.
#
# Every failure is survivable and none of them stop the run: no
# /proc/self/status (macOS, the BSDs, Windows) means the affinity cannot be read
# and there is no taskset there either; an absent or refusing taskset leaves the
# run on the CPUs it had; and a process already pinned to one CPU is left on it,
# so "taskset -c 5 Rscript scale.R" keeps CPU 5 rather than being overruled.
cpus_allowed <- function() {
    path <- "/proc/self/status"
    if (!file.exists(path)) return(NA_character_)
    hit <- grep("^Cpus_allowed_list:", readLines(path, warn = FALSE), value = TRUE)
    if (length(hit) == 0) return(NA_character_)   # a kernel too old to report it
    sub("^Cpus_allowed_list:[[:space:]]*", "", hit[1])
}

pin_to_one_cpu <- function() {
    want <- Sys.getenv("SCALE_CPU", "0")
    if (want == "none" || want == "") {
        cat("SCALE_CPU=none: leaving the scheduler to place the measurements\n")
        return(invisible(NULL))
    }
    if (!grepl("^[0-9]+$", want))
        stop(sprintf("SCALE_CPU must be a CPU number or 'none', not '%s'", want),
             call. = FALSE)

    before <- cpus_allowed()
    if (is.na(before)) {
        cat("cannot read this process's CPU affinity, so the measurements are",
            "not pinned; on a hybrid CPU, start this under the platform's",
            "affinity tool\n")
        return(invisible(NULL))
    }
    # a list with no comma and no dash is one CPU, so someone has already pinned us
    if (!grepl("[,-]", before)) {
        cat(sprintf("already pinned to CPU %s; keeping it\n", before))
        return(invisible(NULL))
    }

    said <- tryCatch(
        suppressWarnings(system2("taskset", c("-pc", want, Sys.getpid()),
                                 stdout = TRUE, stderr = TRUE)),
        error = function(e) conditionMessage(e))
    after <- cpus_allowed()
    if (!is.na(after) && identical(after, want)) {
        cat(sprintf(paste0("pinned to CPU %s (was %s); plot.scaling.pl ",
                           "--measure and scale.py must be on the same one\n"),
                    after, before))
    } else {
        cat(sprintf("could not pin to CPU %s, staying on %s: %s\n",
                    want, before, paste(said, collapse = " ")),
            file = stderr())
    }
    invisible(NULL)
}

pin_to_one_cpu()

results <- list()
too_slow <- character(0)

by_figure <- list()
for (bm in BENCHMARKS) {
    by_figure[[bm$figure]] <- c(by_figure[[bm$figure]], list(bm))
}

for (figure in c("vector", "transform", "io", "frame")) {
    if (is.null(by_figure[[figure]])) next
    for (n in LADDER[[figure]]) {
        n <- as.integer(n)
        todo <- Filter(function(bm) !(bm$name %in% too_slow), by_figure[[figure]])
        if (length(todo) == 0) next
        data <- BUILD[[ BUILDER_FOR[[figure]] ]](n)
        for (bm in todo) {
            # One untimed call.  R resolves methods and grows its vector heap
            # on first use, and this is also where a broken call announces
            # itself before it has been recorded seven times.
            ok <- tryCatch({ bm$body(data); TRUE },
                           error = function(e) { message(sprintf(
                               "%s at n=%d: %s", bm$name, n,
                               conditionMessage(e))); FALSE })
            if (!ok) {
                too_slow <- c(too_slow, bm$name)
                next
            }

            slowest <- 0
            reps <- 0L
            for (run in seq_len(RUNS) - 1L) {
                m <- measure(bm$body, data)
                slowest <- max(slowest, m$seconds)
                reps <- m$reps
                results[[length(results) + 1L]] <- data.frame(
                    figure = figure, `function` = bm$name, call = bm$call,
                    n = n, run = run, seconds = m$seconds,
                    stringsAsFactors = FALSE, check.names = FALSE)
            }
            cat(sprintf("%-9s %-30s n=%-8d %.6f s%s\n", figure, bm$name, n,
                        slowest, if (reps > 1) sprintf(" (x%d)", reps) else ""))
            if (slowest > CAP) too_slow <- c(too_slow, bm$name)
        }
        rm(data)
        invisible(gc(verbose = FALSE))
    }
}

# A curve that ends early is not a curve that was never measured, and the plot
# cannot tell you which one you are looking at.
if (length(too_slow) > 0) {
    cat(sprintf("Stopped early (a run exceeded %g s, or it failed): %s\n",
                CAP, paste(sort(unique(too_slow)), collapse = ", ")))
}

out <- file.path(DIR, sprintf("out.r.%d.tmp", Sys.getpid()))
if (file.exists(out)) unlink(out)

all <- do.call(rbind, results)
write.table(all, "r_scaling.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Done. %d measurements written to r_scaling.tsv\n", nrow(all)))
