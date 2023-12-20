# ================================== CODE METADATA =========================== #
# AUTHORS: Fernando Colchero
# DATE CREATED: 2023-10-01
# DESCRIPTION: Reproducibility code for Sach et al. (2024).
# NOTES: Save and open the  
# ================================ CODE START ================================ #
# ======================== #
# ==== GENERAL SETUP: ====
# ======================== #
# Set working directory:
setwd("path to directory...")

# Source functions:
source("01code/ElephantLongevityFunctions.R")

# Load life tables:
lifeTabs <- read.csv(file = "02data/ElephantLifeTables.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Load number of individuals:
Ninds <- read.csv(file = "02data/ElephantNinds.csv",
                  header = TRUE, stringsAsFactors = FALSE)

# Load Siler mortality parameters:
thetaMat <- read.csv(file = "02data/ElephantSilerParams.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Load life expectancy-lifespan equality:
leMat <- read.csv(file = "02data/ElephantLifeExpLifespEq.csv",
                  header = TRUE, stringsAsFactors = FALSE)

# Plotting colors:
cols <- c("#FD8D3C", "#BD0026", "#377EB8", "#4DAF4A", "#6A3D9A")
names(cols) <- c("ZIMS_1960-1990", "ZIMS_1990-2020", "Amboseli",
                 "Samburu", "Myanmar")

# =================================== #
# ==== PLOT LIFE TABLE SURVIVAL: ====
# =================================== #
# Select species:
sps <- "Loxodonta africana"

# Select sex:
sx <- "Female"

# Subset life tables:
idlt <- which(lifeTabs$Species == sps & lifeTabs$Sex == sx)
spLt <- lifeTabs[idlt, ]

# Find populations:
pops <- unique(spLt$Pop)
npops <- length(pops)

# plot limits:
xlim <- c(0, max(spLt$Ages))
ylim <- c(0, 1)

# Minimum age (e.g., avoid first year):
minAge <- 1

# Prepare empty plot:
par(mfrow = c(1, 1), mar = c(4, 4, 1, 1))
plot(xlim, ylim, col = NA, xlab = "Age", ylab = "Survival")
for (ip in 1:npops) {
  idpop <- which(spLt$Pop == pops[ip])
  ilt <- spLt[idpop, ]
  x <- ilt$Ages[which(ilt$Ages >= minAge)]
  lx <- ilt$lx[which(ilt$Ages >= minAge)]
  lx <- lx / lx[1]
  lines(x, lx, type = 's', col = cols[pops[ip]])
}
legend("topright", legend = pops, col = cols[pops], bty = 'n', lwd = 2)


# =============================================== #
# ==== CALCULATE LIFE EXPECT.-LIFESP. EQUAL. ====
# =============================================== #
# Select species:
sps <- "Loxodonta africana"

# Select sex:
sx <- "Female"

# Select population:
pop <- "ZIMS_1990-2020"

# Subset Siler mortality parameters:
idth <- which(thetaMat$Species == sps & thetaMat$Sex == sx & 
                thetaMat$Pop == pop)
theta <- thetaMat$Mean[idth]
names(theta) <- thetaMat$parameter[idth]

# Age increments (to approximate integral):
dx <- 0.001

# Vector of ages for integration:
x <- seq(0, 200, dx)

# Calculate Survival:
Sx <- CalcSurv(theta = theta, x = x, model = "GO", shape = "bathtub")

# Life expectancy:
ex <- sum(Sx * dx)

# Find ages where Sx > 0:
idn0 <- which(Sx > 0)

# Lifespan equality:
epsx <- -log(-sum(Sx[idn0] * log(Sx[idn0]) * dx) / ex)

# Print results to console:
cat(sprintf("\n============================\nSpecies: %s\nSex: %s\n- Life exp.  = %s\n- Lifesp eq. = %s\n============================\n", sps, sx, round(ex, 2), round(epsx, 3)))

# ================================== #
# ==== KULLBACK-LEIBLER DISCR.: ====
# ================================== #
# Select species:
sps <- "Loxodonta africana"

# Select sex:
sx <- "Female"

# Select population 1:
pop1 <- "ZIMS_1990-2020"

# Select population 2:
pop2 <- "ZIMS_1960-1990"

# Subset leMat:
id1 <- which(leMat$Species == sps & leMat$Sex == sx & leMat$Pop == pop1)
id2 <- which(leMat$Species == sps & leMat$Sex == sx & leMat$Pop == pop2)
le1 <- leMat[id1, ]
le2 <- leMat[id2, ]

# Univariate KL (e.g., life expectancy):
KLuni <- CalcKLc(mean1 = le1$Mean[1], mean2 = le2$Mean[1],
                 sd1 = sqrt(le1$Cov1[1]), sd2 = sqrt(le2$Cov1[1]))

# Calculate bivariate Kullback-Leibler discrepancy:
KLbiv <- CalcBivKL(m1 = le1$Mean, m2 = le2$Mean, 
                cv1 = as.matrix(le1[, c("Cov1", "Cov2")]),
                cv2 = as.matrix(le2[, c("Cov1", "Cov2")]))

