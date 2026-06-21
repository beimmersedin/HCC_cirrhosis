# =============================================================================
#  GRN Network Visualization (igraph/ggraph)
#  GRN.R 결과를 읽어서 네트워크 시각화
# =============================================================================

library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(ggplot2)

# ─── 경로 설정 ───────────────────────────────────────────────────────────────
output_dir <- "/data1/project/yeonu/065_multi_atac/GRN/results_output/"

final_grn      <- read.csv(paste0(output_dir, "Hepatocyte_Final_Weighted_GRN.csv"))
master_summary <- read.csv(paste0(output_dir, "Hepatocyte_Master_Regulators_Summary.csv"))

# ─── Edge / Node Table 생성 ──────────────────────────────────────────────────
edge_table <- final_grn %>%
  group_by(tf_name, gene) %>%
  summarise(weight = max(grn_weight), .groups = "drop") %>%
  rename(source = tf_name, target = gene) %>%
  mutate(interaction = ifelse(weight > 0, "activates", "represses"))

tf_names  <- unique(edge_table$source)
all_nodes <- unique(c(edge_table$source, edge_table$target))

out_degree <- edge_table %>% count(source) %>% rename(node = source, out_degree = n)
in_degree  <- edge_table %>% count(target) %>% rename(node = target, in_degree = n)

node_table <- data.frame(node = all_nodes, stringsAsFactors = FALSE) %>%
  mutate(type = ifelse(node %in% tf_names, "TF", "TargetGene")) %>%
  left_join(out_degree, by = "node") %>%
  left_join(in_degree, by = "node") %>%
  mutate(
    out_degree = ifelse(is.na(out_degree), 0, out_degree),
    in_degree  = ifelse(is.na(in_degree), 0, in_degree),
    total_degree = out_degree + in_degree
  )

top_tfs <- master_summary$tf_name[1:5]


# =============================================================================
# Plot 1: Top 5 Master Regulator Network
# =============================================================================
edges_top5 <- edge_table %>% filter(source %in% top_tfs)
nodes_top5 <- node_table %>% filter(node %in% c(edges_top5$source, edges_top5$target))

graph_top5 <- tbl_graph(nodes = nodes_top5, edges = edges_top5, directed = TRUE) %>%
  mutate(importance = ifelse(type == "TF", out_degree, 1))

set.seed(42)
p1 <- ggraph(graph_top5, layout = 'stress') +
  geom_edge_link(aes(width = weight), alpha = 0.2,
                 arrow = arrow(length = unit(2, 'mm')),
                 end_cap = circle(3, 'mm'), color = "grey") +
  geom_node_point(aes(color = type, size = importance)) +
  geom_node_text(aes(label = ifelse(type == "TF", node, "")),
                 repel = TRUE, size = 5, family = "sans") +
  scale_color_manual(values = c("TF" = "#E41A1C", "TargetGene" = "#377EB8")) +
  scale_size_continuous(range = c(2, 10)) +
  scale_edge_width(range = c(0.2, 1.5)) +
  theme_graph(base_family = "sans") +
  labs(title = "Top 5 Master Regulator Network (Hepatocytes)",
       subtitle = "Validated by TF Activity - Gene Expression Correlation",
       caption = "Edges filtered by Spearman correlation > 0.1")

ggsave(paste0(output_dir, "Hepatocyte_GRN_Top5_Network.png"), plot = p1,
       width = 10, height = 8, dpi = 300)
message("저장 완료: Hepatocyte_GRN_Top5_Network.png")


# =============================================================================
# Plot 2: Full GRN Network
# =============================================================================
graph_all <- tbl_graph(nodes = node_table, edges = edge_table, directed = TRUE) %>%
  mutate(importance = ifelse(type == "TF", out_degree, 1))

top10_tfs <- node_table %>%
  filter(type == "TF") %>%
  arrange(desc(out_degree)) %>%
  head(10) %>%
  pull(node)

set.seed(42)
p2 <- ggraph(graph_all, layout = 'stress') +
  geom_edge_link(aes(width = weight), alpha = 0.1,
                 arrow = arrow(length = unit(1.5, 'mm')),
                 end_cap = circle(2, 'mm'), color = "grey70") +
  geom_node_point(aes(color = type, size = importance), alpha = 0.8) +
  geom_node_text(aes(label = ifelse(node %in% top10_tfs, node, "")),
                 repel = TRUE, size = 4, family = "sans") +
  scale_color_manual(values = c("TF" = "#E41A1C", "TargetGene" = "#377EB8")) +
  scale_size_continuous(range = c(1, 12)) +
  scale_edge_width(range = c(0.1, 1.2)) +
  theme_graph(base_family = "sans") +
  labs(title = "Full GRN Network (Hepatocytes)",
       subtitle = paste0(nrow(node_table), " nodes, ", nrow(edge_table),
                         " edges | Top 10 TF labels shown"),
       caption = "Edges filtered by Spearman correlation > 0.1")

ggsave(paste0(output_dir, "Hepatocyte_GRN_Full_Network.png"), plot = p2,
       width = 14, height = 12, dpi = 300)
message("저장 완료: Hepatocyte_GRN_Full_Network.png")

message("\nDone! All plots saved to: ", output_dir)
