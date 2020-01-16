context("bayz")

test_that("Genotype effect model", {

    n_samples <- 5000

    Genotype_names <- c("A", "B", "C", "D", "E", "F")

    Geno_effects <- c(25, 40, 5, -10, 80, -60)
    names(Geno_effects) <- Genotype_names

    Genotypes <- Genotype_names[runif(n_samples, 1, length(Genotype_names)+1)]

    intercept <- 45
    noise <- rnorm(n_samples, mean=0, sd=30)

    slope <- 0.7

    Response <- intercept + 0.7 * sapply(Genotypes, function(x) Geno_effects[x]) + noise
    names(Response) <- NULL


    fit <- bayz(Response ~ ranf(Genotypes, prior=ichi()), chain=c(2e4,100,10))

    print(summary(fit))

    expect_gte(summary(fit)$Heritability$Heritability[2], 0.5)

})

test_that("Genotype and Environment effects model", {

    n_samples <- 5000

    Genotype_names <- c("AA", "AB", "AC", "AD",
                        "BA", "BB", "BC", "BD",
                        "CA", "CB", "CC", "CD",
                        "DA", "DB", "DC", "DD")

    Geno_effects <- c(200, 150, 70, 30,
                      125, 175, 40, 75,
                      25, -40, -100, 10,
                      80, -25, -80, 125)
    names(Geno_effects) <- Genotype_names

    Genotypes <- Genotype_names[runif(n_samples, 1, length(Genotype_names)+1)]

    Env_1_names <- c("0C", "5C", "10C", "15C", "20C", "25C", "30C")
    Env_1_effects <- c(-100, -40, 10, 160, 60, 20, -30)
    names(Env_1_effects) <- c(Env_1_names)

    Env_1 <- Env_1_names[runif(n_samples, 1, length(Env_1_names)+1)]

    Env_2_names <- c("ph5", "ph6", "ph7", "ph8", "ph9")
    Env_2_effects <- c(-60, 10, 80, 20, -40)
    names(Env_2_effects) <- c(Env_2_names)

    Env_2 <- Env_2_names[runif(n_samples, 1, length(Env_2_names)+1)]

    intercept <- 300
    noise <- rnorm(n_samples, mean=0, sd=40)

    slope <- 0.7

    Response <- intercept + sapply(Genotypes, function(x) Geno_effects[x]) +
    2*sapply(Env_1, function(x) Env_1_effects[x]) +
    2*sapply(Env_2, function(x) Env_2_effects[x])
    names(Response) <- NULL

    fit <- bayz(Response ~ ranf(Genotypes), chain=c(3e4,100,10))

    print(summary(fit))

    expect_gte(summary(fit)$Heritability$Heritability[2], 0.15)
    expect_lte(summary(fit)$Heritability$Heritability[2], 0.35)

    fit <- bayz(Response ~ ranf(Genotypes) + fixf(Env_1) + fixf(Env_2), chain=c(3e4,100,10))

    print(summary(fit))

    expect_gte(summary(fit)$Heritability$Heritability[2], 0.95)

})

test_that("Genotype and Environment interactions effects model", {

    n_samples <- 5000

    Genotype_names <- c("AA", "AB", "AC", "AD",
                        "BA", "BB", "BC", "BD",
                        "CA", "CB", "CC", "CD",
                        "DA", "DB", "DC", "DD")

    Geno_effects <- c(200, 150, 70, 30,
                      125, 175, 40, 75,
                      25, -40, -100, 10,
                      80, -25, -80, 125)
    names(Geno_effects) <- Genotype_names

    Genotypes <- Genotype_names[runif(n_samples, 1, length(Genotype_names)+1)]

    Env_1_names <- c("0C", "5C", "10C", "15C", "20C", "25C", "30C")
    Env_1_effects <- c(-100, -40, 10, 160, 60, 20, -30)
    names(Env_1_effects) <- c(Env_1_names)

    Env_1 <- Env_1_names[runif(n_samples, 1, length(Env_1_names)+1)]

    Env_2_names <- c("ph5", "ph6", "ph7", "ph8", "ph9")
    Env_2_effects <- c(-60, 10, 80, 20, -40)
    names(Env_2_effects) <- c(Env_2_names)

    Env_2 <- Env_2_names[runif(n_samples, 1, length(Env_2_names)+1)]

    intercept <- 300
    noise <- rnorm(n_samples, mean=0, sd=40)

    slope <- 0.7

    Response <- intercept + sapply(Genotypes, function(x) Geno_effects[x]) +
    2*sapply(Env_1, function(x) Env_1_effects[x]) +
    2*sapply(Env_2, function(x) Env_2_effects[x]) +
    0.03*(sapply(Env_1, function(x) Env_1_effects[x])*sapply(Env_2, function(x) Env_2_effects[x]))
    names(Response) <- NULL

    fit <- bayz(Response ~ ranf(Genotypes), chain=c(3e4,100,10))

    print(summary(fit))

    expect_gte(summary(fit)$Heritability$Heritability[2], 0.05)
    expect_lte(summary(fit)$Heritability$Heritability[2], 0.25)

    fit <- bayz(Response ~ ranf(Genotypes) + fixf(Env_1) + fixf(Env_2), chain=c(3e4,100,10))

    print(summary(fit))

    expect_gte(summary(fit)$Heritability$Heritability[2], 0.4)
    expect_lte(summary(fit)$Heritability$Heritability[2], 0.6)

})
