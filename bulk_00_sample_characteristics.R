setwd("")

library(dplyr)


df <- read.csv("")

check <- df %>%
  group_by(GSE_id) %>%
  summarise(
    across(
      c(group, age, sex, pmi, batch),
      ~ n_distinct(.x, na.rm = FALSE)
    ),
    .groups = "drop"
  )

if (any(as.matrix(check[, -1]) > 1)) {
  stop("Inconsistent metadata within the same donor.")
}

df <- df %>% distinct(GSE_id, .keep_all = TRUE)
df$group <- factor(df$group, levels = c("CTRL", "BD"))

table(df$group)
# CTRL   BD 
# 39   32 


test_numeric <- function(data, variable) {
  
  x <- data[data$group == "CTRL", variable]
  y <- data[data$group == "BD", variable]
  
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  p <- if (length(x) > 0 && length(y) > 0) {
    wilcox.test(x, y, exact = FALSE)$p.value
  } else {
    NA_real_
  }
  
  summary_fun <- function(z) {
    if (length(z) == 0) return(NA_character_)
    
    sprintf(
      "%.2f [%.2f–%.2f]",
      median(z),
      quantile(z, 0.25),
      quantile(z, 0.75)
    )
  }
  
  data.frame(
    Variable = variable,
    Test = "Wilcoxon rank-sum",
    CTRL = summary_fun(x),
    BD = summary_fun(y),
    P_value = p
  )
}

test_categorical <- function(data, variable) {
  
  sub <- data[
    !is.na(data[[variable]]) &
      !is.na(data$group),
  ]
  
  tab <- table(
    sub$group,
    sub[[variable]]
  )
  
  p <- if (nrow(tab) == 2 && ncol(tab) >= 2 &&
           all(rowSums(tab) > 0)) {
    tryCatch(
      fisher.test(tab)$p.value,
      error = function(e) {
        fisher.test(
          tab,
          simulate.p.value = TRUE,
          B = 100000
        )$p.value
      }
    )
  } else {
    NA_real_
  }
  
  summary_fun <- function(group_name) {
    
    if (!(group_name %in% rownames(tab))) {
      return(NA_character_)
    }
    
    counts <- tab[group_name, ]
    total <- sum(counts)
    
    if (total == 0) return(NA_character_)
    
    paste(
      sprintf(
        "%s: %d (%.1f%%)",
        names(counts),
        as.integer(counts),
        100 * counts / total
      ),
      collapse = "; "
    )
  }
  
  data.frame(
    Variable = variable,
    Test = "Fisher's exact",
    CTRL = summary_fun("CTRL"),
    BD = summary_fun("BD"),
    P_value = p
  )
}

results <- bind_rows(
  test_numeric(df, "age"),
  test_categorical(df, "sex"),
  test_numeric(df, "pmi"),
  test_categorical(df, "batch")
)

results$FDR <- p.adjust(results$P_value, method = "BH")

results$Significance <- ifelse(
  is.na(results$FDR),
  NA_character_,
  ifelse(results$FDR < 0.05, "Yes", "No")
)

print(results)

write.csv(results, "", row.names = FALSE)
