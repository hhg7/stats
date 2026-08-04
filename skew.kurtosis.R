# Set random seed for reproducibility
set.seed(42)
n <- 10000

# 1. BASELINE (Normal)
norm_data <- rnorm(n, mean = 0, sd = 1)

# 2. SKEWNESS
right_skew_data <- rlnorm(n, meanlog = 0, sdlog = 0.6)  # Log-Normal
left_skew_data  <- rbeta(n, shape1 = 8, shape2 = 2)     # Beta (requires shifting/scaling to overlay perfectly, but shows the shape)

# 3. KURTOSIS
heavy_tail_data <- rt(n, df = 3)                        # Student's t (Leptokurtic)
light_tail_data <- runif(n, min = -3, max = 3)          # Uniform (Platykurtic)

# --- PLOTTING ---
par(mfrow=c(1,2)) # Set up a 1x2 plot area

# Plot 1: Skewness
plot(density(norm_data), main="Visualizing Skewness", lwd=2, ylim=c(0, 0.8), xlim=c(-3, 6))
lines(density(right_skew_data), col="red", lwd=2)
lines(density(left_skew_data), col = "blue", lwd=2)
legend("topright", legend=c("Normal (Baseline)", "Positive Skew (Log-Normal)"), 
       col=c("black", "red"), lwd=2)

# Plot 2: Kurtosis
plot(density(norm_data), main="Visualizing Kurtosis", lwd=2, ylim=c(0, 0.5), xlim=c(-6, 6))
lines(density(heavy_tail_data), col="blue", lwd=2)
lines(density(light_tail_data), col="green", lwd=2)
legend("topright", legend=c("Normal (Mesokurtic)", "Heavy Tails (t-dist)", "Light Tails (Uniform)"), 
       col=c("black", "blue", "green"), lwd=2)
