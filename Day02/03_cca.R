set.seed(42)

# Simulate 100 (genes) x 150 (cells) matrix with Gaussian noise
d1 <- matrix(rnorm(100 * 150, mean = 0, sd = 1), nrow = 100, ncol = 150)

# Create d2 by adding 10 to all entries
d2 <- d1 + 10

create_heatmap <- function(mat, title, filename=NULL) {
    if (!is.null(filename)){
        png(filename, width = 800, height = 600, res = 100)
    }

    colors <- colorRampPalette(c("blue", "white", "red"))(100)

    breaks <- seq(min(mat), max(mat), length.out = 101)

    image(t(mat)[, nrow(mat):1],
          col = colors,
          breaks = breaks,
          main = title,
          xlab = "Column",
          ylab = "Row",
          axes = FALSE)

    axis(1, at = seq(0, 1, length.out = 5),
         labels = round(seq(1, ncol(mat), length.out = 5)))
    axis(2, at = seq(0, 1, length.out = 5),
         labels = round(seq(nrow(mat), 1, length.out = 5)))
    if (!is.null(filename)){

    dev.off()
    }
}

# Create visualizations of the matrices
create_heatmap(d1, "Matrix d1 (Gaussian noise, mean=0, sd=1)")

create_heatmap(d1, "Matrix d1 (Gaussian Noise, mean=0, sd=1)",
               "matrix_d1.png")
create_heatmap(d2, "Matrix d2 (d1 + 10)")

create_heatmap(d2, "Matrix d2 (d1 + 10)",
               "matrix_d2.png")

# Perform CCA
cca_result <- cancor(d1, d2)

# Extract canonical correlations
canonical_cors <- cca_result$cor

plot(1:length(canonical_cors), canonical_cors,
     type = "b",
     col = "steelblue",
     lwd = 2,
     pch = 19,
     cex = 1.5,
     main = "Canonical correlations between d1 and d2",
     xlab = "Canonical cariate pair",
     ylab = "Canonical correlation",
     ylim = c(0, 1))
#grid()
#dev.off()

# The cancor function returns coefficient matrices of dimension min(p,q) x min(p,q)
# where p and q are the number of columns in d1 and d2

# Check dimensions of coefficient matrices
cat("\nDimensions of xcoef:", dim(cca_result$xcoef), "\n")
cat("Dimensions of ycoef:", dim(cca_result$ycoef), "\n")

# Create a combined plot with 3 panels showing the relationship
par(mfrow = c(1, 3))

# For visualization, we'll plot samples from d1 vs d2 for first 3 dimensions
for (i in 1:min(3, min(ncol(d1), 100))) {
    # Use columns from the original matrices
    plot(d1[, i], d2[, i],
         pch = 19,
         col = rgb(0.2, 0.4, 0.8, 0.5),
         main = paste0("Dimension ", i, " Comparison\n",
                       "d1 vs d2 (shifted by +10)"),
         xlab = paste0("d1 Column ", i),
         ylab = paste0("d2 Column ", i))

    # Add regression line
    abline(lm(d2[, i] ~ d1[, i]), col = "red", lwd = 2)
    grid()
}

#dev.off()

# All canonical correlations are essentially 1.0 because d2 = d1 + 10


# Now let's play with different d2

d2 <- matrix(0, nrow = 100, ncol = 150)

# Columns 1-50: Strong correlation (0.9)
for (i in 1:50) {
    d2[, i] <- 0.9 * d1[, i] + 0.1 * rnorm(100, mean = 0.1, sd = 2)
}

# Columns 51-100: Moderate correlation (0.5)
for (i in 51:100) {
    d2[, i] <- 0.5 * d1[, i] + 0.5 * rnorm(100, mean = 0.1, sd = 2)
}

# Columns 101-150: Weak correlation (mostly noise)
for (i in 101:150) {
    d2[, i] <- 0.01 * d1[, i] + 1 * rnorm(100, mean = 0.1, sd = 2)
}

n_rows <- 200
n_cols <- 50

# Simulate n_rows x n_cols matrix with Gaussian noise
d1 <- matrix(rnorm(n_rows * n_cols, mean = 0, sd = 1), nrow = n_rows, ncol = n_cols)

# Create d2 with varying correlation structure
# First 15 columns: strong correlation (0.9)
# Next 15 columns: moderate correlation (0.6)
# Next 10 columns: weak correlation (0.3)
# Last 10 columns: independent noise

d2 <- matrix(rnorm(n_rows * n_cols, mean = 0, sd = 1), nrow = n_rows, ncol = n_cols)

# Apply correlations
for (i in 1:15) {
    d2[, i] <- 0.9 * d1[, i] + sqrt(1 - 0.9^2) * d2[, i]
}

for (i in 16:30) {
    d2[, i] <- 0.6 * d1[, i] + sqrt(1 - 0.6^2) * d2[, i]
}

for (i in 31:40) {
    d2[, i] <- 0.3 * d1[, i] + sqrt(1 - 0.3^2) * d2[, i]
}

create_heatmap(d1, "Matrix d1")
create_heatmap(d2, "Matrix d2 new")


cca_result <- cancor(d1, d2)

# Extract canonical correlations
canonical_cors <- cca_result$cor

plot(1:length(canonical_cors), canonical_cors,
     type = "b",
     col = "steelblue",
     lwd = 2,
     pch = 19,
     cex = 1.5,
     main = "Canonical correlations between d1 and d2",
     xlab = "Canonical cariate pair",
     ylab = "Canonical correlation",
     ylim = c(0, 1))
