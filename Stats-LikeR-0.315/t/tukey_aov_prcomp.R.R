# Regenerate the frozen table in t/tukey_aov_prcomp.R.t:
#
#     Rscript t/tukey_aov_prcomp.R.R
#
# and paste the output into that file's data section.  Written against
# R 4.6.1 (2026-06-24).  The test reads only the frozen literals; installing or
# testing Stats::LikeR needs no R.
#
# Provenance -- all of it R's own material:
#
#   ptukey / qtukey  <- stats::ptukey, stats::qtukey.  LikeR.xs carries a
#       faithful port of R's src/nmath/{ptukey,qtukey}.c (Copenhaver & Holland
#       1988), so the right test is a dense grid against R itself rather than a
#       handful of textbook table entries.  637 points: nmeans in
#       {2,3,4,5,6,10,20} x df in {3,5,10,27,60,120,1e6}, at seven quantiles
#       and six probabilities.  t/tukey.t checks a few of these to 1e-3, which
#       understates the port by eleven orders of magnitude -- see the .t file.
#
#   aov / TukeyHSD   <- stats::aov, stats::TukeyHSD
#       src/library/stats/man/PlantGrowth.Rd and its use in TukeyHSD's own
#       examples: weight ~ group on PlantGrowth, 3 groups of 10.
#       datasets::chickwts: weight ~ feed, 6 groups of unequal size, which is
#       what exercises the unbalanced studentized-range interval.
#
#   anova            <- stats::anova on an lm, sequential Type I
#       mtcars, mpg ~ wt + hp + disp -- three terms, so term order matters and
#       a Type-I/Type-III mix-up shows up immediately.
#
#   vif              <- 1/(1 - R^2) of each predictor on the others, which is
#       what car::vif computes; done here with base lm() so the generator needs
#       no contributed package.
#
#   prcomp           <- stats::prcomp
#       src/library/stats/man/prcomp.Rd uses USArrests; both scale.=TRUE and
#       the unscaled default are emitted.
#
#   scale            <- base::scale, including R's center=FALSE divisor (the
#       root mean square about zero, not the standard deviation).
#
# Output is one row per case: label <TAB> comma-separated values.

options(digits = 17, warn = -1)

f <- function(k, v) cat(sprintf("%s\t%s\n", k,
	paste(sapply(v, function(z)
		if (is.null(z) || length(z) == 0 || is.na(z)) "NaN" else sprintf("%.17g", z)),
	      collapse = ",")))

## ---- ptukey / qtukey grid ----
for (nm in c(2, 3, 4, 5, 6, 10, 20)) {
	for (df in c(3, 5, 10, 27, 60, 120, 1e6)) {
		for (q in c(0.5, 1, 2, 3, 3.5, 5, 8))
			f(sprintf("ptukey|%g|%g|%g", q, nm, df), ptukey(q, nm, df))
		for (p in c(0.5, 0.75, 0.9, 0.95, 0.99, 0.999))
			f(sprintf("qtukey|%g|%g|%g", p, nm, df), qtukey(p, nm, df))
	}
}
# the upper tail, and the identity qtukey(p, upper) == qtukey(1-p, lower)
f("ptukey.upper|3|3|27",   ptukey(3, 3, 27, lower.tail = FALSE))
f("qtukey.upper|0.05|3|27", qtukey(0.05, 3, 27, lower.tail = FALSE))
# df -> Inf is sqrt(2) * qnorm(1 - alpha/2) for nmeans = 2, which is the one
# closed form the studentized range has
f("qtukey.inf2", qtukey(0.95, 2, 1e9))
f("qnorm.975.sqrt2", sqrt(2) * qnorm(0.975))

## ---- aov + TukeyHSD on PlantGrowth ----
data(PlantGrowth)
a <- aov(weight ~ group, data = PlantGrowth)
s <- summary(a)[[1]]
f("aov.pg", c(s[1, "Df"], s[1, "Sum Sq"], s[1, "Mean Sq"], s[1, "F value"],
              s[1, "Pr(>F)"], s[2, "Df"], s[2, "Sum Sq"], s[2, "Mean Sq"]))
th <- TukeyHSD(a)$group
for (i in seq_len(nrow(th))) f(paste0("tukey.pg|", rownames(th)[i]), th[i, ])
for (cl in c(0.9, 0.99)) {
	th2 <- TukeyHSD(a, conf.level = cl)$group
	for (i in seq_len(nrow(th2)))
		f(sprintf("tukey.pg.cl%g|%s", cl * 100, rownames(th2)[i]), th2[i, ])
}
cat(sprintf("# PlantGrowth weight: %s\n", paste(PlantGrowth$weight, collapse = ",")))
cat(sprintf("# PlantGrowth group:  %s\n", paste(PlantGrowth$group, collapse = ",")))

## ---- aov + TukeyHSD on chickwts (6 unbalanced groups) ----
data(chickwts)
a <- aov(weight ~ feed, data = chickwts)
s <- summary(a)[[1]]
f("aov.cw", c(s[1, "Df"], s[1, "Sum Sq"], s[1, "Mean Sq"], s[1, "F value"],
              s[1, "Pr(>F)"], s[2, "Df"], s[2, "Sum Sq"], s[2, "Mean Sq"]))
th <- TukeyHSD(a)$feed
for (i in seq_len(nrow(th))) f(paste0("tukey.cw|", rownames(th)[i]), th[i, ])
for (lv in levels(chickwts$feed))
	cat(sprintf("# chickwts %s = %s\n", lv,
	            paste(chickwts$weight[chickwts$feed == lv], collapse = ",")))

## ---- anova: sequential Type I on mtcars ----
data(mtcars)
m  <- lm(mpg ~ wt + hp + disp, data = mtcars)
an <- anova(m)
for (i in seq_len(nrow(an))) {
	fv <- if (is.na(an[i, "F value"])) NA else an[i, "F value"]
	pv <- if (is.na(an[i, "Pr(>F)"])) NA else an[i, "Pr(>F)"]
	f(paste0("anova.mt|", rownames(an)[i]),
	  c(an[i, "Df"], an[i, "Sum Sq"], an[i, "Mean Sq"], fv, pv))
}
f("lm.mt.coef", coef(m))
f("lm.mt.fit",  c(summary(m)$r.squared, summary(m)$adj.r.squared, summary(m)$sigma))
# vif: 1/(1 - R^2) from regressing each predictor on the other two
f("vif.mt", c(1 / (1 - summary(lm(wt   ~ hp + disp, mtcars))$r.squared),
              1 / (1 - summary(lm(hp   ~ wt + disp, mtcars))$r.squared),
              1 / (1 - summary(lm(disp ~ wt + hp,   mtcars))$r.squared)))
for (v in c("mpg", "wt", "hp", "disp"))
	cat(sprintf("# mtcars %s = %s\n", v, paste(mtcars[[v]], collapse = ",")))

## ---- prcomp on USArrests (prcomp.Rd's own example) ----
data(USArrests)
for (sc in c(TRUE, FALSE)) {
	tag <- if (sc) "prcomp.sc" else "prcomp.ns"
	p <- prcomp(USArrests, scale. = sc)
	f(paste0(tag, ".sdev"), p$sdev)
	for (j in 1:4) f(sprintf("%s.rot%d", tag, j), p$rotation[, j])
	if (sc) { f(paste0(tag, ".center"), p$center); f(paste0(tag, ".scale"), p$scale) }
	else      f(paste0(tag, ".center"), p$center)
}
cat("# USArrests column order: Murder,Assault,UrbanPop,Rape\n")
for (i in seq_len(nrow(USArrests)))
	cat(sprintf("# usa %s\n", paste(as.numeric(USArrests[i, ]), collapse = ",")))

## ---- scale ----
f("scale.1to5",     as.numeric(scale(1:5)))
f("scale.nocenter", as.numeric(scale(1:5, center = FALSE)))
f("scale.noscale",  as.numeric(scale(1:5, scale = FALSE)))
f("scale.neither",  as.numeric(scale(1:5, center = FALSE, scale = FALSE)))
f("scale.ctr3",     as.numeric(scale(1:5, center = 3)))
f("scale.ctr3.sc2", as.numeric(scale(1:5, center = 3, scale = 2)))
mm <- matrix(c(1, 3, 5, 2, 4, 6), ncol = 2)          # columns (1,3,5) and (2,4,6)
f("scale.mat",    as.numeric(scale(mm)))
f("scale.mat.nc", as.numeric(scale(mm, center = FALSE)))
f("scale.mat.ns", as.numeric(scale(mm, scale = FALSE)))
