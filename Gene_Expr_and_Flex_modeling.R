library(dplyr)
library(tidyr)
library(nlme)
library(purrr)
library(ggplot2)
library(future)
library(furrr)
library(progressr)
library(edgeR)


####Read input data and format####
all <- read.csv("acclim_data_all.csv") #dataset contains all individual names and experimental information
counts <- read.csv(file = "/merged_counts_all_111422.csv") #this is a merged featureCounts table including all taxa
rownames(counts) <- counts$X
counts$X <- NULL
counts <- na.omit(counts)

all <- all[all$Mouse_ID %in% colnames(counts),]

counts$mean = rowMeans(counts) #rowMeans takes means of each row
keep_counts = subset(counts, mean >= 10) #filter out genes that have <10 reads on average across individuals
dim(keep_counts)
keep_counts$mean = NULL #clean up dataset
dim(keep_counts)

treatment <- all$Treatment #this is my first factor for the GLM model
species <- all$Strain #this is my second factor
table <- data.frame(treatment, species)
group <- factor(paste(table$treatment, table$species, sep="_"))
y <- DGEList(counts=keep_counts, group=group) # make a DGE list
dim(y)

y <- calcNormFactors(y) # normalize
plotMDS(y) #exploration plots
counts.norm <- cpm(y, log=TRUE, prior.count=2, normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset
#write.csv(counts.norm, file = "11Mar25_cpm_counts.csv")

counts.norm <- as.data.frame(counts.norm)

cpm_t <- t(counts.norm)
#format input df
cpm_meta <- all[all$Mouse_ID %in% rownames(cpm_t),]
match(cpm_meta$Mouse_ID, rownames(cpm_t))

cpm_expr_meta <- cpm_t
cpm_expr_meta$Strain <- cpm_meta$Strain
cpm_expr_meta$Treatment <- cpm_meta$Treatment
cpm_expr_meta$Elevation <- cpm_meta$Elevation
cpm_expr_meta$Family <- cpm_meta$Family

#reformat df so we can scale up to multiple species
expr_long <- cpm_expr_meta %>%
  pivot_longer(
    cols = -c(Treatment, Strain, Elevation),
    names_to  = "gene",
    values_to = "expression"
  )

strain_elev <- tibble::tibble(
  Strain = c("PO", "LS", "LN", "ML", "AZ", "ME"),
  Elevation = c(0, 430, 430, 1500, 2700, 4350)
)

#build function for gls model
fit_gls_by_strain <- function(df) {
  mod <- gls(
    expression ~ Treatment,
    data = df,
    na.action = na.omit
  )
  
  tibble(
    estimate  = coef(mod)["Treatment"],
    std_error = summary(mod)$tTable["Treatment", "Std.Error"],
    t_value   = summary(mod)$tTable["Treatment", "t-value"],
    p_value   = summary(mod)$tTable["Treatment", "p-value"]
  )
}

#test one gene
A4galt_res <- expr_long %>%
  filter(gene == "A4galt") %>%
  group_by(Strain) %>%
  group_modify(~ fit_gls_by_strain(.x)) %>%
  ungroup()

A4galt_res <- A4galt_res %>%
  left_join(strain_elev, by = "Strain")

A4galt_pgls_model <- nlme::gls(
  estimate ~ Elevation,
  data = A4galt_res,
  correlation = ape::corPagel(1, CAE_cytb_tree, form = ~ Strain)
)

A4galt_slope_summary <- summary(A4galt_pgls_model)

A4galt_pgls_estimate <- coef(A4galt_pgls_model)["Elevation"]
A4galt_pgls_pval     <- summary(A4galt_pgls_model)$tTable["Elevation", "p-value"]

fit_pgls_on_slopes <- function(df) {
  mod <- nlme::gls(
    estimate ~ Elevation,
    data = df,
    correlation = ape::corPagel(1, CAE_cytb_tree, form = ~ Strain)
  )
  
  tibble(
    slope_elev   = coef(mod)["Elevation"],
    p_value      = summary(mod)$tTable["Elevation", "p-value"],
    lambda       = mod$modelStruct$corStruct[1]
  )
}

A4galt_pgls_res <- A4galt_res %>%
  fit_pgls_on_slopes()

fit_data <- data.frame(
  elev = hct_scaled_df$elev,
  fitted_values = predict(slope_model)
)

#scale up gls (Expr ~ Treatment | Family)function for all genes
fit_gene_treatment_gls <- function(data, gene) {
  
  gls_data <- data %>%
    dplyr::select(Strain, Treatment, Family, !!gene) %>%
    dplyr::rename(expr = !!gene) %>%
    dplyr::mutate(
    # expr = arm::rescale(expr),
    TreatVal = as.numeric(as.character(Treatment))
     )
  
  gls_data %>%
    dplyr::group_by(Strain) %>%
    dplyr::summarise(
      gls_model = list(
        nlme::gls(
          expr ~ TreatVal,
          #weights     = nlme::varIdent(form = ~ 1 | Treatment),
          #correlation = nlme::corCompSymm(form = ~ 1 | Family),
          na.action   = na.omit,
          data        = dplyr::pick(dplyr::everything())
        )
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gene      = gene,
      estimate  = purrr::map_dbl(gls_model, ~ coef(.x)["TreatVal"]),
      std_error = purrr::map_dbl(gls_model, ~ summary(.x)$tTable["TreatVal", "Std.Error"]),
      t_value   = purrr::map_dbl(gls_model, ~ summary(.x)$tTable["TreatVal", "t-value"]),
      p_value   = purrr::map_dbl(gls_model, ~ summary(.x)$tTable["TreatVal", "p-value"])
    ) %>%
    dplyr::select(gene, Strain, estimate, std_error, t_value, p_value)
}

#filter genes with at least one taxon with a significant repsonse
gene_passes_filter <- function(gls_df, alpha = 0.05) {
  gls_df %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(any_sig = any(p_value < alpha, na.rm = TRUE)) %>%
    dplyr::filter(any_sig) %>%
    dplyr::pull(gene)
}

#PGLS per treatment-sig gene
fit_gene_pgls <- function(gls_df) {
  
  pgls_data <- gls_df %>%
    dplyr::left_join(strain_elev, by = "Strain") %>%
    dplyr::mutate(Strain = factor(Strain, levels = CAE_cytb_tree$tip.label))
  
  model <- nlme::gls(
    estimate ~ Elevation,
    data = pgls_data,
    correlation = ape::corPagel(1, CAE_cytb_tree, form = ~ Strain)
  )
  
  tibble::tibble(
    gene        = unique(pgls_data$gene),
    slope       = coef(model)["Elevation"],
    p_value     = summary(model)$tTable["Elevation", "p-value"],
    lambda      = model$modelStruct$corStruct[1]
  )
}

####run GLS across all genes####
genes <- setdiff(
  colnames(cpm_expr_meta),
  c("Treatment", "Strain", "Elevation", "Family")
)

treatment_gls_all <- purrr::map_dfr(
  genes,
  ~ fit_gene_treatment_gls(cpm_expr_meta, .x)
)

#running into problems with false convergence, 
#where sample size is too low, or not enough within-treatment variance
#try a safety wrapper, where genes that pass get shuttled to next step and genes
#that don't are reported
safe_fit_gene_treatment_gls <- purrr::safely(
  fit_gene_treatment_gls,
  otherwise = NULL
)

treatment_gls_list <- purrr::map(
  genes,
  ~ safe_fit_gene_treatment_gls(cpm_expr_meta, .x)
)

# Successful fits, formatted for easy viewing
treatment_gls_all <- treatment_gls_list %>%
  purrr::keep(~ !is.null(.x$result)) %>%
  purrr::map_dfr("result")
length(unique(treatment_gls_all$gene))

#and this is the pgls_input format, with the gls estimates retained for each gene
treatment_gls_list_filt <- purrr::keep(
  treatment_gls_list,
  ~ !is.null(.x$result) &&
    any(.x$result$p_value < 0.05, na.rm = TRUE)
)

#clean this up a bit
has_sig_treatment <- function(x, alpha = 0.05) {
  !is.null(x$result) &&
    any(x$result$p_value < alpha, na.rm = TRUE)
}

treatment_gls_list_filt <- purrr::keep(
  treatment_gls_list,
  has_sig_treatment
)


# Failed genes
failed_genes <- treatment_gls_list %>%
  purrr::keep(~ !is.null(.x$error)) %>%
  purrr::map_chr(~ .x$error$message)
length(failed_genes)

####run sig filter####
genes_to_test <- gene_passes_filter(treatment_gls_all)

filtered_gls <- treatment_gls_all %>%
  dplyr::filter(gene %in% genes_to_test)
length(unique(filtered_gls$gene))

####run pgls on filtered genes####
#First, build func for indiv gene and add Elev to list
fit_gene_pgls_abs <- function(gene_res) {
  
  df <- gene_res$result %>%
    dplyr::left_join(strain_elev, by = "Strain") %>%
    dplyr::mutate(
      abs_estimate = abs(estimate)
    )
  
  #Filter because PGLS requires variation
  if (nrow(df) < 3 || dplyr::n_distinct(df$Elevation) < 3) {
    stop("Insufficient elevation variation for PGLS")
  }
  
  model <- nlme::gls(
    abs_estimate ~ Elevation,
    data = df,
    correlation = ape::corPagel(1, CAE_cytb_tree, form = ~ Strain)
  )
  
  smry <- summary(model)$tTable
  
  tibble::tibble(
    gene        = unique(df$gene),
    slope       = smry["Elevation", "Value"],
    std_error   = smry["Elevation", "Std.Error"],
    t_value     = smry["Elevation", "t-value"],
    p_value     = smry["Elevation", "p-value"]
  )
}




#similar convegrence error as in GLS above, safety wrap it
safe_fit_gene_pgls_abs <- purrr::safely(
  fit_gene_pgls_abs,
  otherwise = NULL
)

pgls_results_abs <- purrr::map(
  treatment_gls_list_filt,
  ~ safe_fit_gene_pgls_abs(.x)
)

pgls_results_keep_abs <- purrr::keep(
    pgls_results_abs,
    ~ !is.null(.x$result)
  )
length(pgls_results_keep_abs)
#this is now a list of genes where the pgls successfully converged (1732 of 1897)

pgls_failed_genes <- pgls_results_abs %>%
  purrr::keep(~ !is.null(.x$error)) %>%
  purrr::imap_dfr(
    ~ tibble::tibble(
      gene  = .y,
      error = .x$error$message
    )
  )

#bind the list of successful PGLS models into one DF
pgls_df_abs <- purrr::map_dfr(
  pgls_results_keep_abs,
  "result"
)


####FDR pval correction####
pgls_df_abs <- pgls_df_abs %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

pgls_sig_abs <- pgls_df_abs %>%
  dplyr::filter(p_adj < 0.05)
dim(pgls_sig_abs)

assimilation_candidates <- pgls_sig_abs %>%
  filter(slope < 0) %>%
  arrange(slope)

assimilation_candidates_clean <- fill_na_from_column("./PGLS_Gene_Expr/9Feb26_reduced_flex_candidates_orthology.csv", "ortholog_name", "initial_alias")
assimilation_candidates_clean <- deduplicate_by_column(assimilation_candidates_clean, "initial_alias")
write.csv(assimilation_candidates_clean, file = "./PGLS_Gene_Expr/9Feb26_reduced_flex_candidates_orthology.csv", row.names = FALSE)

####PGLS Slope Visualization####
#use earlier tibble list and join with Elev for input df
treatment_gls_with_elev <- purrr::map(
  treatment_gls_list_filt,
  ~ .x$result %>%
    dplyr::left_join(strain_elev, by = "Strain") %>%
    dplyr::rename(slope = estimate)
)

treatment_gls_with_elev <- treatment_gls_with_elev %>%
  rlang::set_names(
    purrr::map_chr(., ~ unique(.x$gene))
  )


input_df <- treatment_gls_with_elev[["St3gal2"]]

#ggplot
#fix small color matching issue
#strain_colors <- setNames(color_frame$color, color_frame$species)

plot_gene_flexibility_pgls_abs <- function(input_df, color_frame,
                                       ylims = NULL) {
  
  # --- safety checks ---
  stopifnot(all(c("Strain", "Elevation", "slope") %in% names(input_df)))
  stopifnot(dplyr::n_distinct(input_df$Elevation) > 1)
  
  # --- abs value of slope --- 
  input_df <- input_df %>%
    dplyr::mutate(abs_slope = abs(slope))
  
  # --- fit PGLS ---
  pgls_model <- nlme::gls(
    abs_slope ~ Elevation,
    data = input_df,
    correlation = ape::corPagel(1, CAE_cytb_tree, form = ~ Strain)
  )
  
  pgls_summary <- summary(pgls_model)
  
  slope_estimate <- coef(pgls_model)[["Elevation"]]
  p_val <- pgls_summary$tTable["Elevation", "p-value"]
  
  # --- fitted line ---
  fit_data <- data.frame(
    Elevation = input_df$Elevation,
    fitted_values = predict(pgls_model)
  )
  
  # --- plotting ---
  p <- ggplot(input_df, aes(x = Elevation, y = abs_slope, color = Strain)) +
    geom_point(size = 6) +
    geom_line(
      data = fit_data,
      aes(x = Elevation, y = fitted_values),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    scale_color_manual(values = strain_colors) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      y = "Flexibility (abs of treatment slope)",
      x = "Native Elevation (m a.s.l.)"
    ) +
    annotate(
      "text",
      x = 3000,
      y = max(input_df$abs_slope, na.rm = TRUE),
      label = paste0(
        "Slope = ", signif(slope_estimate, 3),
        "\np = ", signif(p_val, 3)
      ),
      hjust = 0,
      size = 4
    )
  
  if (!is.null(ylims)) {
    p <- p + coord_cartesian(ylim = ylims)
  }
  
  list(
    plot = p,
    model = pgls_model,
    summary = pgls_summary
  )
}

####plot abs flex ~ elevation with a single candidate gene####
gene_df <- treatment_gls_with_elev[["Srf"]]
gene_flex_pgls <- plot_gene_flexibility_pgls_abs(
  input_df = gene_df,
  color_frame = color_frame,
  #ylims = c(-2e-4, 2e-4)   # optional, but recommended for comparability
)

gene_flex_pgls$plot
gene_flex_pgls$model
gene_flex_pgls$summary

####Split by Treatment, Run Summary Statistics, and fit PGLS####
#split into a named list
treatment_dfs <- split(cpm_expr_meta, cpm_expr_meta$Treatment)
#convert to long format
expr_long_by_treatment <- purrr::map(
  treatment_dfs,
  ~ .x %>%
    pivot_longer(
      cols = -c(Strain, Elevation, Treatment),
      names_to = "gene",
      values_to = "expression"
    )
)
#gate to gene being expressed > 2 per strain, calc summary stats
summarise_gene_expression <- function(df, min_n = 3) {
  df %>%
    group_by(gene, Strain, Elevation) %>%
    summarise(
      n = sum(!is.na(expression)),
      mean_expression = mean(expression, na.rm = TRUE),
      se_expression   = sd(expression, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) %>%
    filter(n >= min_n)
}
#apply across all treatments
gene_means_by_treatment <- purrr::map(
  expr_long_by_treatment,
  summarise_gene_expression
)
#due to PGLS requiring variation, we need to ensure each gene is represented by all 6 taxa,
#but this can be adjusted
filter_genes_with_enough_strains <- function(df, min_strains = 6) {
  df %>%
    group_by(gene) %>%
    filter(dplyr::n_distinct(Strain) >= min_strains) %>%
    ungroup()
}
gene_means_by_treatment_filt <- purrr::map(
  gene_means_by_treatment,
  filter_genes_with_enough_strains
)

#build pgls function
fit_gene_pgls <- function(df, gene_id, tree) {
  
  gene_df <- df %>% filter(gene == gene_id)
  
  # extra safety gate
  if (nrow(gene_df) < 4) return(NULL)
  
  model <- nlme::gls(
    mean_expression ~ Elevation,
    data = gene_df,
    correlation = ape::corPagel(1, tree, form = ~ Strain)
  )
  
  tibble(
    gene      = gene_id,
    slope     = coef(model)["Elevation"],
    std_error = summary(model)$tTable["Elevation", "Std.Error"],
    t_value   = summary(model)$tTable["Elevation", "t-value"],
    p_value   = summary(model)$tTable["Elevation", "p-value"]
  )
}

safe_fit_gene_pgls <- purrr::safely(
  fit_gene_pgls,
  otherwise = NULL
)

run_pgls_for_treatment <- function(df, tree) {
  
  genes <- unique(df$gene)
  
  purrr::map(
    genes,
    ~ safe_fit_gene_pgls(df, .x, tree)
  ) %>%
    purrr::keep(~ !is.null(.x$result)) %>%
    purrr::map_dfr("result")
}

pgls_results_by_treatment <- purrr::map(
  gene_means_by_treatment,
  run_pgls_for_treatment,
  tree = CAE_cytb_tree
)

pgls_results_by_treatment <- purrr::map(
  pgls_results_by_treatment,
  ~ .x %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
)

pgls_results_by_treatment$`1000`
pgls_results_by_treatment$`3000`
pgls_results_by_treatment$`4500`

converg_fail_1000 <- setdiff(unique(gene_means_by_treatment[["1000"]]$gene), 
        pgls_results_by_treatment[["1000"]]$gene)

converg_fail_3000 <- setdiff(unique(gene_means_by_treatment[["3000"]]$gene), 
                             pgls_results_by_treatment[["3000"]]$gene)

converg_fail_4500 <- setdiff(unique(gene_means_by_treatment[["4500"]]$gene), 
                             pgls_results_by_treatment[["4500"]]$gene)

# Ensure uniqueness within each vector (important if duplicates are possible)
converg_fail_1000 <- unique(converg_fail_1000)
converg_fail_3000 <- unique(converg_fail_3000)
converg_fail_4500 <- unique(converg_fail_4500)

# Unique genes per treatment (present in one, absent in others)
unique_1000 <- setdiff(converg_fail_1000, union(converg_fail_3000, converg_fail_4500))
unique_3000 <- setdiff(converg_fail_3000, union(converg_fail_1000, converg_fail_4500))
unique_4500 <- setdiff(converg_fail_4500, union(converg_fail_1000, converg_fail_3000))

# Percent unique relative to each treatment
length(unique_1000) / length(converg_fail_1000)
length(unique_3000) / length(converg_fail_3000)
length(unique_4500) / length(converg_fail_4500)

# (Optional) Percent relative to all genes observed across treatments
all_genes <- unique(c(converg_fail_1000, converg_fail_3000, converg_fail_4500))

pct_all_1000 <- count_1000 / length(all_genes)
pct_all_3000 <- count_3000 / length(all_genes)
pct_all_4500 <- count_4500 / length(all_genes)

# Combine into a clean summary table
summary_df <- data.frame(
  treatment = c("1000", "3000", "4500"),
  unique_gene_count = c(count_1000, count_3000, count_4500),
  pct_within_treatment = c(pct_1000, pct_3000, pct_4500),
  pct_of_all_genes = c(pct_all_1000, pct_all_3000, pct_all_4500)
)

summary_df

# Count significant genes and save to a separate list
sig_genes_by_treatment <- purrr::map(
  pgls_results_by_treatment,
  ~ .x %>%
    filter(p_adj < 0.05) %>%
    pull(gene)
)
# Quick summary
purrr::map_int(sig_genes_by_treatment, length)

#plot mean expression ~ elevation
plot_candidate_gene <- function(df, gene_id, tree, color_frame) {
  
  # Subset gene
  gene_df <- df %>% filter(gene == gene_id)
  
  # Fit PGLS
  pgls_model <- nlme::gls(
    mean_expression ~ Elevation,
    data = gene_df,
    correlation = ape::corPagel(1, tree, form = ~ Strain)
  )
  
  # Prepare fitted line
  fit_data <- gene_df %>%
    mutate(fitted = predict(pgls_model))
  
  # Extract slope and p-value for annotation
  slope_estimate <- coef(pgls_model)["Elevation"]
  p_val <- summary(pgls_model)$tTable["Elevation", "p-value"]
  
  # Plot
  ggplot(gene_df, aes(x = Elevation, y = mean_expression, color = Strain)) +
    geom_point(size = 5) +
    geom_line(data = fit_data, aes(x = Elevation, y = fitted), 
              color = "black", linewidth = 0.5, linetype = "dashed") +
    scale_color_manual(values = strain_colors) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      y = "Mean Expression",
      x = "Native Elevation (m a.s.l.)"
    ) +
    annotate(
      "text", 
      x = min(gene_df$Elevation) + 200, 
      y = max(gene_df$mean_expression, na.rm = TRUE),
      label = paste0("Slope = ", round(slope_estimate, 8), 
                     "\nP = ", round(p_val, 8)),
      hjust = 0, size = 4, color = "black"
    )
}
