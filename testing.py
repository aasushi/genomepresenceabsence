#223
gtdb_ids = list(con.index)
working_dir = os.getcwd()

ko_cols = [c for c in con.columns if c.startswith('K')]
con[ko_cols].to_csv(os.path.join(working_dir, "ko_trait_matrix.csv"))

r_script = f"""
library(ape)
library(phytools)
tree <- read.tree("bac120_r214.tree")
taxa_to_keep <- c({','.join(f'"{g}"' for g in gtdb_ids)})
taxa_to_keep <- intersect(tree$tip.label, taxa_to_keep)
pruned_tree <- keep.tip(tree, taxa_to_keep)

traits <- read.csv("ko_trait_matrix.csv", row.names = 1)
shared <- intersect(pruned_tree$tip.label, rownames(traits))
pruned_tree <- keep.tip(pruned_tree, shared)
traits <- traits[shared, ]
traits <- traits[, apply(traits, 2, var) > 0]
ppca <- phyl.pca(pruned_tree, traits, method = "BM", mode = "corr")
write.csv(ppca$L, "ppca_loadings.csv")
write.csv(ppca$S, "ppca_scores.csv")
eig_df <- data.frame(
    PC = rownames(ppca$Eval),
    eigenvalue = diag(ppca$Eval),
    variance_explained = diag(ppca$Eval) / sum(diag(ppca$Eval)) * 100)
write.csv(eig_df, "ppca_eigenvalues.csv", row.names = FALSE)
pdf("ppca_biplot.pdf", width=10, height=10)
biplot(ppca)
dev.off()
cat("pPCA complete\\n")
"""

with open("prune_tree.R", "w") as f:
    f.write(r_script)

subprocess.run(["Rscript", "prune_tree.R"], check=True)
print("Tree pruned and saved")

ppca_scores = pd.read_csv("ppca_scores.csv", index_col=0)
ppca_loadings = pd.read_csv("ppca_loadings.csv", index_col=0)
ppca_eigenvalues = pd.read_csv("ppca_eigenvalues.csv")



os.makedirs("itol_annotations", exist_ok=True)

ko_colors = {
    'K00370': '#e41a1c',
    'K02567': '#377eb8',
    'K00368': '#4daf4a',
    'K15864': '#984ea3',
    'K04561': '#ff7f00',
    'K00376': '#a65628',
    'K00411': '#f781bf',
    'K00412': '#999999',
    'K00413': '#ffff33',}

ko_labels = {
    'K00370': 'NarG',
    'K02567': 'Nap',
    'K00368': 'NirK',
    'K15864': 'NirS',
    'K04561': 'cNor',
    'K00376': 'NosZ',
    'K00411': 'bc1-subunit1',
    'K00412': 'bc1-subunit2',
    'K00413': 'bc1-subunit3',}


for ko, label in ko_labels.items():
    lines = [
        "DATASET_BINARY",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t{label}",
        f"COLOR\t{ko_colors[ko]}",
        "FIELD_LABELS\tpresence",
        "FIELD_SHAPES\t1",
        "DATA",]

    for genome_id in con.index:
        val = con.loc[genome_id, ko]
        lines.append(f"{genome_id}\t{val}")

    with open(f"itol_annotations/{ko}_{label}.txt", "w") as f:
        f.write("\n".join(lines))

#for assembly status
status_colors = {
    'complete/high quality': '#2ecc71',
    'near complete':         '#f39c12',
    'medium quality draft':  '#e67e22',
    'fragmented':            '#e74c3c'}

lines = [
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tAssembly Status",
    "COLOR\t#000000",
    "LEGEND_TITLE\tAssembly Status",
    "LEGEND_SHAPES\t1\t1\t1\t1",
    "LEGEND_COLORS\t"+"\t".join(status_colors.values()),
    "LEGEND_LABELS\t"+"\t".join(status_colors.keys()),
    "DATA"]

for genome_id in con.index:
    status = con.loc[genome_id, 'assembly_status']
    color = status_colors.get(status, '#cccccc')
    lines.append(f"{genome_id}\t{color}\t{status}")

with open("itol_annotations/assembly_status.txt", "w") as f:
    f.write("\n".join(lines))

#for denitrification steps
step_colors = {0: '#ffffff', 1: '#fee8c8', 2: '#fdbb84', 3: '#e34a33', 4: '#b30000'}

lines = [
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tDenitrification Steps",
    "COLOR\t#000000",
    "LEGEND_TITLE\tSteps Encoded",
    "LEGEND_SHAPES\t1\t1\t1\t1\t1",
    "LEGEND_COLORS\t" + "\t".join(step_colors.values()),
    "LEGEND_LABELS\t0 steps\t1 step\t2 steps\t3 steps\t4 steps",
    "DATA"]

for genome_id in con.index:
    steps = 0
    if (con.loc[genome_id, step_one] == 1).any():
        steps += 1
    if (con.loc[genome_id, step_two] == 1).any():
        steps += 1
    if con.loc[genome_id, step_three] == 1:
        steps += 1
    if con.loc[genome_id, step_four] == 1:
        steps += 1
    color = step_colors[steps]
    lines.append(f"{genome_id}\t{color}\t{steps} steps")

with open("itol_annotations/denitrification_steps.txt", "w") as f:
    f.write("\n".join(lines))

print("iTOL annotation files saved")


