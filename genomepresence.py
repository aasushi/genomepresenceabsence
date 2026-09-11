import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os

csvs = ["K00368_Annotree_Hits.csv","K00370_Annotree_Hits.csv","K00376_Annotree_Hits.csv",
        "K02567_Annotree_Hits.csv","K02575_Annotree_Hits.csv","K03385_Annotree_Hits.csv",
        "K04561_Annotree_Hits.csv","K07234_Annotree_Hits.csv","K07673_Annotree_Hits.csv",
        "K10535_Annotree_Hits.csv","K12266_Annotree_Hits.csv","K13771_Annotree_Hits.csv",
        "K15864_Annotree_Hits.csv","K21563_Annotree_Hits.csv","K21564_Annotree_Hits.csv",
        "K00411_Annotree_Hits.csv","K00412_Annotree_Hits.csv","K00413_Annotree_Hits.csv"]

df = list()

for csv in csvs:
    try:
        d = pd.read_csv(csv)
        #split name at _ and take everything from 0th index to then
        KO_number = csv.split('_')[0]
        #drop unneeded rows
        dropped = ['geneId', 'sequence', 'tax','SearchId']
        d = d.drop(dropped, axis = 1)
        #make the KO_number the column labels
        d[KO_number] = 1
        #make gtdbId the row labels
        d = d[['gtdbId', KO_number]].set_index('gtdbId')
        dfinal = d.groupby('gtdbId')[KO_number].max().to_frame()
        #put dataframes into new list
        df.append(dfinal)
        #print("CSV file imported successfully!")
    except:
        print(f"Error: The file '{csv}' was not found.")

#combine dataframes and fill any unfilled values with 0
if df:
    con = pd.concat(df,axis=1)
    con = con.fillna(0)
    con = con.astype(int)

#add new columns for contamination and completeness
con['contamination'] = None
con['completeness'] = None
con['assembly_status'] = None

contam = "bac120_metadata_r214.tsv"
chunksize = 10**2
for chunk in pd.read_csv(contam, sep='\t', chunksize = chunksize):
    for index, row in chunk.iterrows():
        genome_id = row['gtdb_genome_representative']
        if genome_id in con.index:
            #locate matching gtdb row, find contamination/compeleteness column, and add bac120 info
            con.loc[genome_id, 'contamination'] = row['checkm_contamination']
            con.loc[genome_id, 'completeness'] = row['checkm_completeness']

#creating the assembly status column with conditions
con['assembly_status'] = 'near complete'
con.loc[(con['completeness'] >= 95) & (con['contamination'] <= 5), 'assembly_status'] = 'complete/high quality'
con.loc[(con['completeness'] >= 90) & (con['completeness'] < 95) & (con['contamination'] <= 5), 'assembly_status'] = 'near complete'
con.loc[(con['completeness'] >= 70) & (con['completeness'] < 90) & (con['contamination'] <= 10), 'assembly_status'] = 'medium quality draft'
con.loc[(con['completeness'] < 70) | (con['contamination'] > 10), 'assembly_status'] = 'fragmented'

#status = con['assembly_status'].value_counts()
#plot = status.plot.pie(figsize =(5,5), autopct='%1.1f%%')

step_one = ['K00370','K02567']
step_two = ['K15864','K00368']
step_three = 'K04561'
step_four = 'K00376'

total_genomes = len(con)
#print(f"Total genomes pre-filter: {total_genomes}")

all_steps = []
for value in con.index:
    step_count = 0
    if (con.loc[value, step_one] == 1).any():
        step_count = step_count + 1
    if (con.loc[value, step_two] == 1).any():
        step_count = step_count + 1
    if con.loc[value, step_three] == 1:
        step_count = step_count + 1
    if con.loc[value, step_four] == 1:
        step_count = step_count + 1
    all_steps.append(step_count)

initial_step_counts = pd.Series(all_steps).value_counts().sort_index()
#print(initial_step_counts)

# List of quality levels to exclude
exclude_levels = [
    ['fragmented'],
    ['fragmented', 'medium quality draft'],
    ['fragmented', 'medium quality draft', 'near complete']]


for excluded in exclude_levels:
    steps_per_genome = list()
    configurations = list()
    for value in con.index:
        if con.loc[value, 'assembly_status'] not in excluded:
            step_count = 0
            config = ''

            if (con.loc[value, step_one] == 1).any():
                step_count = step_count + 1
                config = config + '1'
            else:
                config = config + '0'

            if (con.loc[value, step_two] == 1).any():
                step_count = step_count + 1
                config = config + '1'
            else:
                config = config + '0'

            if con.loc[value, step_three] == 1:
                step_count = step_count + 1
                config = config + '1'
            else:
                config = config + '0'

            if con.loc[value, step_four] == 1:
                step_count = step_count + 1
                config = config + '1'
            else:
                config = config + '0'

            configurations.append(config)
            steps_per_genome.append(step_count)

    n_genomes = len(configurations)
    reduction = ((total_genomes - n_genomes)/total_genomes)*100

    #print(f"Excluding: {excluded}")
    #print(f"Genomes remaining: {n_genomes}")
    #print(f"Genomes excluded: {total_genomes - n_genomes}")
    #print(f"Reduction: {reduction:.1f}%")

    step_counts = pd.Series(steps_per_genome).value_counts().sort_index()
    #print(step_counts)
    #print(f"Step count reductions from initial:")
    for step in initial_step_counts.index:
        if step in step_counts.index:
            initial = initial_step_counts[step]
            current = step_counts[step]
            step_reduction = ((initial - current)/initial)*100
            #print(f"{step} steps:({step_reduction:.1f}% reduction)")

    plot = step_counts.plot.pie(figsize=(5, 5), autopct='%1.1f%%')
    plt.title(f"Number of Denitrification Steps Encoded per Genome excluding {excluded}")
    plt.ylabel('')
    #plt.show()

    config_count = pd.Series(configurations).value_counts()
    #print(f"Configuration Counts excluding {excluded}")

    plot = config_count.plot.pie(figsize=(8, 8), autopct='%1.1f%%')
    plt.title(f'Genome Configurations Excluding {excluded}')
    plt.ylabel('')
   # plt.show()

bc1genes = ["K00411","K00412","K00413"]
downstream_kos = ['K02567', 'K00368', 'K15864', 'K04561', 'K00376', 'K00370']

filter_levels = [
    'No filter',
    'Excl. fragmented',
    'Excl. frag + med',
    'Only complete']
filter_excluded = [
    [],
    ['fragmented'],
    ['fragmented', 'medium quality draft'],
    ['fragmented', 'medium quality draft', 'near complete']]

reductase_names = ['Nap', 'NirK', 'NirS', 'cNor', 'NosZ', 'NarG']

proportions = []
rows = []

for excluded, filter_label in zip(filter_excluded, filter_levels):
    if excluded:
        sub = con[~con['assembly_status'].isin(excluded)]
    else:
        sub = con

    #print(f"{filter_label}: {len(sub)} genomes")

    for name, ko in zip(reductase_names, downstream_kos):
        has_reductase = sub[sub[ko] == 1]
        n = len(has_reductase)

        # of those, how many also have bc1
        has_bc1 = (has_reductase[bc1genes] == 1).all(axis=1).sum()
        pct = (has_bc1/n)*100

        #print(f"{name}: {n} genomes with reductase, {has_bc1} also have bc1 ({pct:.1f}%)")
        rows.append({'filter': filter_label, 'reductase': name, 'n_with_reductase': n, 'n_with_bc1': int(has_bc1),
                     'pct_with_bc1': round(pct, 1)})

    has_reductase = sub[(sub[downstream_kos] == 1).any(axis=1)]
    has_bc1 = (has_reductase[bc1genes] == 1).all(axis=1).sum()

    n = len(has_reductase)
    pct = (has_bc1/n)*100
    proportions.append(round(pct, 1))
    #print(f"  Overall: {n} genomes with any reductase, {has_bc1} also have bc1 ({pct:.1f}%)")

results_df = pd.DataFrame(rows)
pivot = results_df.pivot(index='reductase', columns='filter', values='pct_with_bc1')
pivot = pivot[filter_levels]

ax = pivot.plot.bar(figsize=(10, 6), width=0.7)
plt.title('% of Genomes with Downstream Reductase that also Encode bc1')
plt.xlabel('Reductase')
plt.ylabel('% also encoding bc1')
plt.legend(title='Quality Filter', loc='upper left')
plt.savefig('bc1_cooccurrence.png')
#plt.show()

for seed in [1, 2, 3]:
    con_filtered = con[con['assembly_status'] == 'complete/high quality']
    con_filtered = con_filtered.sample(n=3000, random_state=seed)
    print(f"Seed {seed} - Genomes used for pPCA: {len(con_filtered)}")

    gtdb_ids = list(con_filtered.index)
    working_dir = os.getcwd()

    ko_cols = [c for c in con_filtered.columns if c.startswith('K')]
    con_filtered[ko_cols].to_csv(os.path.join(working_dir, "ko_trait_matrix.csv"))

    r_script = f"""
library(ape)
library(phytools)
setwd("{working_dir}")

tree <- read.tree("bac120_r214.tree")
taxa_to_keep <- c({','.join(f'"{g}"' for g in gtdb_ids)})
taxa_to_keep <- intersect(tree$tip.label, taxa_to_keep)
pruned_tree <- keep.tip(tree, taxa_to_keep)
write.tree(pruned_tree, "pruned_bac120_seed{seed}.tree")

traits <- read.csv("ko_trait_matrix.csv", row.names = 1)
shared <- intersect(pruned_tree$tip.label, rownames(traits))
pruned_tree <- keep.tip(pruned_tree, shared)
traits <- traits[shared, ]

name_map <- c(K00370="NarG", K02567="Nap", K00368="NirK", K15864="NirS",
              K04561="cNor", K00376="NosZ", K00411="bc1-sub1", K00412="bc1-sub2",
              K00413="bc1-sub3", K02575="NarK", K03385="NrfA", K07234="NosL",
              K07673="NorD", K10535="Hao", K12266="NosD", K13771="NarK2",
              K21563="NosF", K21564="NosY")
colnames(traits) <- ifelse(colnames(traits) %in% names(name_map),
                           name_map[colnames(traits)],
                           colnames(traits))

traits <- traits[, apply(traits, 2, var) > 0]
ppca <- phyl.pca(pruned_tree, traits, method = "BM", mode = "corr")
write.csv(ppca$L, "ppca_loadings_seed{seed}.csv")
write.csv(ppca$S, "ppca_scores_seed{seed}.csv")
eig_df <- data.frame(
    PC = rownames(ppca$Eval),
    eigenvalue = diag(ppca$Eval),
    variance_explained = diag(ppca$Eval) / sum(diag(ppca$Eval)) * 100)
write.csv(eig_df, "ppca_eigenvalues_seed{seed}.csv", row.names = FALSE)
pdf("ppca_biplot_seed{seed}.pdf", width=14, height=14)
biplot(ppca,
       cex = c(0.0001, 0.8),
       col = c("grey70", "red"),
       arrow.len = 0.05)
dev.off()
cat("pPCA seed {seed} complete\\n")
"""

    with open(f"ppca_seed{seed}.R", "w") as f:
        f.write(r_script)

    env = os.environ.copy()
    env["R_MAX_VSIZE"] = "32Gb"
    subprocess.run(["Rscript", f"ppca_seed{seed}.R"], check=True, env=env)
    print(f"Seed {seed} complete!")

    ppca_scores = pd.read_csv(f"ppca_scores_seed{seed}.csv", index_col=0)
    ppca_loadings = pd.read_csv(f"ppca_loadings_seed{seed}.csv", index_col=0)
    ppca_eigenvalues = pd.read_csv(f"ppca_eigenvalues_seed{seed}.csv")
    ppca_scores = ppca_scores.join(con_filtered[['assembly_status']])

    # Scree plot
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(ppca_eigenvalues['PC'], ppca_eigenvalues['variance_explained'], color='steelblue')
    ax.set_xlabel('PC')
    ax.set_ylabel('Variance Explained (%)')
    ax.set_title(f'Bar Plot - Seed {seed}')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'bar_plot_seed{seed}.png', dpi=150)
    plt.close()
    print(f"Bar plot seed {seed} saved!")

    # Step count scatter plot
    step_counts = []
    for genome_id in con_filtered.index:
        steps = 0
        if (con_filtered.loc[genome_id, step_one] == 1).any():
            steps += 1
        if (con_filtered.loc[genome_id, step_two] == 1).any():
            steps += 1
        if con_filtered.loc[genome_id, step_three] == 1:
            steps += 1
        if con_filtered.loc[genome_id, step_four] == 1:
            steps += 1
        step_counts.append(steps)

    con_filtered = con_filtered.copy()
    con_filtered['step_count'] = step_counts
    ppca_scores = ppca_scores.join(con_filtered[['step_count']])

    step_colors = {0: '#ffffff', 1: '#fee8c8', 2: '#fdbb84', 3: '#e34a33', 4: '#b30000'}
    fig, ax = plt.subplots(figsize=(8, 8))
    for steps, group in ppca_scores.groupby('step_count'):
        ax.scatter(group['PC1'], group['PC2'],
                   label=f'{steps} steps',
                   color=step_colors[steps],
                   edgecolors='grey',
                   linewidths=0.3,
                   alpha=0.8,
                   s=30)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.legend(title='Denitrification Steps')
    plt.title(f'Phylogenetic PCA - Seed {seed} (colored by denitrification steps)')
    plt.savefig(f'ppca_scatter_steps_seed{seed}.png', dpi=150)
    plt.close()
    print(f"Step count scatter plot seed {seed} saved!")

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

# for assembly status
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

# for denitrification steps
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