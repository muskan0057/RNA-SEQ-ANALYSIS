
library(GSA)
library(gprofiler2)
library(writexl)
library(openxlsx)


#g1 UP in InfVsUI_NT but not in InfVsUI_Plin2
#g2 UP in InfVsUI_Plin2 but not in InfVsUI_NT
#g3 DOWN in InfVsUI_NT but not in InfVsUI_Plin2
#g4 DOWN in InfVsUI_Plin2 but not in InfVsUI_NT
g1 = read.xlsx("UP_NTsh.xlsx")
g2 = read.xlsx("UP_Plin2sh.xlsx")
g3 = read.xlsx("DN_NTsh.xlsx")
g4 = read.xlsx("DN_Plin2sh.xlsx")
g5 = read.xlsx("UP_Common.xlsx")
g6 = read.xlsx("DN_Common.xlsx")
g11 = as.vector(g1)
g22 = as.vector(g2)
g33 = as.vector(g3)
g44 = as.vector(g4)
g55 = as.vector(g5)
g66 = as.vector(g6)

gmt7 = upload_GMT_file("c7.all.v2023.2.Hs.symbols.gmt")
gmt2 = upload_GMT_file("c2.cp.v2023.2.Hs.symbols.gmt")
gmt5 = upload_GMT_file("c5.all.v2023.2.Hs.symbols.gmt")
gmt3 = upload_GMT_file("c3.all.v2023.2.Hs.symbols.gmt")
gmtH = upload_GMT_file("h.all.v2023.2.Hs.symbols.gmt")

g1_gmt3 = gost(
  query = g11,
  organism = gmt3[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g1_gmtH = gost(
  query = g11,
  organism = gmtH[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g2_gmt3 = gost(
  query = g22,
  organism = gmt3[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g2_gmtH = gost(
  query = g22,
  organism = gmtH[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g3_gmt3 = gost(
  query = g33,
  organism = gmt3[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g3_gmtH = gost(
  query = g33,
  organism = gmtH[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g4_gmt3 = gost(
  query = g44,
  organism = gmt3[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g4_gmtH = gost(
  query = g44,
  organism = gmtH[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g5_gmt3 = gost(
  query = g55,
  organism = gmt3[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)


g5_gmtH = gost(
  query = g55,
  organism = gmtH[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g6_gmt3 = gost(
  query = g66,
  organism = gmt3[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g6_gmtH = gost(
  query = g66,
  organism = gmtH[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)

g5_gmt2 = gost(
  query = g55,
  organism = gmt2[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE)

g5_gmt5 = gost(
  query = g55,
  organism = gmt5[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE)

g5_gmt7 = gost(
  query = g55,
  organism = gmt7[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE
)


g6_gmt2 = gost(
  query = g66,
  organism = gmt2[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE)

g6_gmt5 = gost(
  query = g66,
  organism = gmt5[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE)

g6_gmt7 = gost(
  query = g66,
  organism = gmt7[[1]],
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE,
  highlight = FALSE)


library(writexl)
setwd("C:/Users/acer/Documents/phd/NTsh_Plin2sh/results")
knitr::opts_knit$set(root.dir = "C:/Users/acer/Documents/phd/NTsh_Plin2sh/results")

writexl::write_xlsx(g1_gmt3$result, "g1_gmt3_results.xlsx")
writexl::write_xlsx(g5_gmt3$result, "g5_gmt3_results.xlsx")


writexl::write_xlsx(g1_gmtH$result, "g1_gmtH_results.xlsx")
writexl::write_xlsx(g2_gmtH$result, "g2_gmtH_results.xlsx")
writexl::write_xlsx(g5_gmtH$result, "g5_gmtH_results.xlsx")
writexl::write_xlsx(g6_gmtH$result, "g6_gmtH_results.xlsx")

writexl::write_xlsx(g5_gmt2$result, "g5_gmt2_results.xlsx")
writexl::write_xlsx(g5_gmt5$result, "g5_gmt5_results.xlsx")
writexl::write_xlsx(g5_gmt7$result, "g5_gmt7_results.xlsx")
writexl::write_xlsx(g6_gmt2$result, "g6_gmt2_results.xlsx")
writexl::write_xlsx(g6_gmt5$result, "g6_gmt5_results.xlsx")
writexl::write_xlsx(g6_gmt7$result, "g6_gmt7_results.xlsx")


gp1_gmt2 = gostplot(g1_gmt2, capped = TRUE, interactive = TRUE)
gp1_gmt7 = gostplot(g1_gmt7, capped = TRUE, interactive = TRUE)
gp2_gmt2 = gostplot(g2_gmt2, capped = TRUE, interactive = TRUE)
gp2_gmt7 = gostplot(g2_gmt7, capped = TRUE, interactive = TRUE)
gp2_gmt5 = gostplot(g2_gmt5, capped = TRUE, interactive = TRUE)
gp3_gmt7 = gostplot(g3_gmt7, capped = TRUE, interactive = TRUE)
gp4_gmt7 = gostplot(g4_gmt7, capped = TRUE, interactive = TRUE)
gostplot(g1_gmt2, capped = TRUE, interactive = TRUE)
gostplot(g1_gmt7, capped = TRUE, interactive = TRUE)
gostplot(g2_gmt2, capped = TRUE, interactive = TRUE)
gostplot(g2_gmt5, capped = TRUE, interactive = TRUE)
gostplot(g2_gmt7, capped = TRUE, interactive = TRUE)
gostplot(g3_gmt7, capped = TRUE, interactive = TRUE)
gostplot(g4_gmt7, capped = TRUE, interactive = TRUE)


gp1_gmt2 = gostplot(g1_gmt2, capped = TRUE, interactive = TRUE)
gp1_gmt7 = gostplot(g1_gmt7, capped = TRUE, interactive = TRUE)
gp2_gmt2 = gostplot(g2_gmt2, capped = TRUE, interactive = TRUE)
gp2_gmt7 = gostplot(g2_gmt7, capped = TRUE, interactive = TRUE)
gp2_gmt5 = gostplot(g2_gmt5, capped = TRUE, interactive = TRUE)
gp3_gmt7 = gostplot(g3_gmt7, capped = TRUE, interactive = TRUE)
gp4_gmt7 = gostplot(g4_gmt7, capped = TRUE, interactive = TRUE)
gostplot(g1_gmt2, capped = TRUE, interactive = TRUE)
gostplot(g1_gmt7, capped = TRUE, interactive = TRUE)
gostplot(g2_gmt2, capped = TRUE, interactive = TRUE)
gostplot(g2_gmt5, capped = TRUE, interactive = TRUE)
gostplot(g2_gmt7, capped = TRUE, interactive = TRUE)
gostplot(g3_gmt7, capped = TRUE, interactive = TRUE)
gostplot(g4_gmt7, capped = TRUE, interactive = TRUE)
