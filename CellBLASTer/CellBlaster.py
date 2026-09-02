import argparse
import hashlib
import os
import re
from pathlib import Path
from urllib.parse import quote
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
import random
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN
from sklearn.neighbors import LocalOutlierFactor
from sklearn.model_selection import train_test_split
from scipy.stats import fisher_exact
from scipy.sparse import issparse
from collections import defaultdict, deque
from numba import config as numba_config
from numba import njit, prange, set_num_threads, get_num_threads

DEFAULT_NONCODING_RNA_KEYWORDS = [
    # Long non-coding RNA
    "LNC",
    "lncRNA",
    "lincRNA",
    "long_non_coding",
    "long-non-coding",
    "long non-coding",
    "long_noncoding",
    "long-noncoding",
    "long noncoding",
    "ncRNA",
    "noncoding_RNA",
    "non-coding_RNA",
    "noncoding RNA",
    "non-coding RNA",
    # MicroRNA
    "miRNA",
    "microRNA",
    "MIR",
    # Small nuclear/nucleolar RNA
    "snoRNA",
    "small_nucleolar_RNA",
    "small nucleolar RNA",
    "snRNA",
    "small_nuclear_RNA",
    "small nuclear RNA",
    "scaRNA",
    "scRNA",
    "snlRNA",
    # Transfer and ribosomal RNA
    "tRNA",
    "transfer_RNA",
    "transfer RNA",
    "rRNA",
    "ribosomal_RNA",
    "ribosomal RNA",
    "Mt_tRNA",
    "Mt_rRNA",
    # Other small and regulatory RNA
    "misc_RNA",
    "miscRNA",
    "siRNA",
    "phasiRNA",
    "tasiRNA",
    "ta-siRNA",
    "antisense_RNA",
    "antisense RNA",
    "antisense",
    "sense_intronic",
    "sense_overlapping",
    "3prime_overlapping_ncRNA",
    # Other non-coding RNA classes
    "RNase_P_RNA",
    "RNase_MRP_RNA",
    "SRP_RNA",
    "signal_recognition_particle_RNA",
    "telomerase_RNA",
    "guide_RNA",
    "vault_RNA",
    "Y_RNA",
    "ribozyme",
    "circRNA",
    "processed_transcript",
]

DATABASE_CONFIG = {
    ("Root", "Monocot"): {
        "record_id": "22233404",
        "orthogroups": "Root_Monocot_Orthogroups.txt",
        "celltype": "Root_celltype.txt",
    },
    ("Root", "Dicot"): {
        "record_id": "21915931",
        "orthogroups": "Root_Dicot_Orthogroups.txt",
        "celltype": "Root_celltype.txt",
    },
    ("Leaf", "Monocot"): {
        "record_id": "22232205",
        "orthogroups": "Leaf_monocot_Orthogroups.txt",
        "celltype": "Leaf_celltype.txt",
    },
    ("Leaf", "Dicot"): {
        "record_id": "22232205",
        "orthogroups": "Leaf_dicot_Orthogroups.txt",
        "celltype": "Leaf_celltype.txt",
    },
    ("Flower", "Monocot"): {
        "record_id": "22232358",
        "orthogroups": "Flower_monocot_Orthogroups.txt",
        "celltype": "Flower_celltype.txt",
    },
    ("Flower", "Dicot"): {
        "record_id": "22232358",
        "orthogroups": "Flower_Dicot_Orthogroups.txt",
        "celltype": "Flower_celltype.txt",
    },
}


class CellBlaster:
    def __init__(
        self,
        database_type,
        symbols,
        output_path,
        query,
        query_symbol,
        filter_keywords=None,
        reference_adata=None,
        reference_symbol=None,
        organ="Root",
        n_jobs=None,
    ):
        if (organ, database_type) not in DATABASE_CONFIG:
            raise ValueError("organ must be Root/Leaf/Flower, and database_type must be Dicot/Monocot")
        self.database_type = database_type
        # Retain legacy attribute names to avoid breaking code that uses the old interface.
        self.dabase_type = database_type
        self.organ = organ
        self.symbols = list(dict.fromkeys(symbols))
        self.output_path = str(Path(output_path).expanduser().resolve())
        self.query = str(Path(query).expanduser().resolve())
        safe_name = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
        invalid_symbols = [symbol for symbol in self.symbols if not safe_name.fullmatch(symbol)]
        if invalid_symbols:
            raise ValueError("Dataset symbols may contain only letters, numbers, periods, underscores, and hyphens")
        if not safe_name.fullmatch(query_symbol):
            raise ValueError("query_symbol may contain only letters, numbers, periods, underscores, and hyphens")
        if (reference_adata is None) != (reference_symbol is None):
            raise ValueError("reference adata and reference symbol must be provided simultaneously")
        if reference_symbol is not None:
            if not safe_name.fullmatch(reference_symbol):
                raise ValueError("reference_symbol may contain only letters, numbers, periods, underscores, and hyphens")
            if reference_symbol in self.symbols:
                raise ValueError("The user reference symbol cannot be duplicated with the downloaded reference symbol")
        self.query_symbol = query_symbol
        self.reference_adata = (
            str(Path(reference_adata).expanduser().resolve())
            if reference_adata is not None
            else None
        )
        self.reference_symbol = reference_symbol
        self.filter_keywords = (
            DEFAULT_NONCODING_RNA_KEYWORDS.copy()
            if filter_keywords is None
            else list(filter_keywords)
        )
        self.parent_graph = defaultdict(set)
        available_cpus = os.cpu_count() or 1
        requested_cpus = n_jobs or int(os.environ.get("BLAST_N_JOBS", "30"))
        self.n_jobs = max(1, min(requested_cpus, available_cpus, numba_config.NUMBA_NUM_THREADS))
        set_num_threads(self.n_jobs)
        print(f"Number of CPU cores used for parallel computation: {get_num_threads()}")

    @staticmethod
    def _checksum_matches(path, checksum):
        if not checksum:
            return True
        algorithm, expected = checksum.split(":", 1)
        digest = hashlib.new(algorithm)
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest().lower() == expected.lower()

    def download_single_data(self, url, dest_dir, filename, checksum=None):
        """Download a single file using a temporary file and checksum to avoid incomplete data."""
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        target_path = dest_dir / filename
        if target_path.is_file():
            print(f"File already exists; skipping download: {filename}")
            return str(target_path)
        # if (target_path.exists() and target_path.stat().st_size > 0 and self._checksum_matches(target_path, checksum)):
        #     print(f"File already exists and passed checksum verification; skipping: {filename}")
        #     return str(target_path)
        temporary_path = target_path.with_suffix(target_path.suffix + ".part")
        print(f"Downloading {filename} ...")
        try:
            with requests.get(url, stream=True, timeout=(30, 300)) as response:
                response.raise_for_status()
                with open(temporary_path, "wb") as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            if not self._checksum_matches(temporary_path, checksum):
                raise RuntimeError(f"Downloaded file failed checksum verification: {filename}")
            os.replace(temporary_path, target_path)
        except Exception as exc:
            temporary_path.unlink(missing_ok=True)
            raise RuntimeError(f"Failed to download {filename}: {exc}") from exc
        print(f"Download completed: {target_path}")
        return str(target_path)

    @staticmethod
    def _find_record_file(record_files, candidates):
        for candidate in candidates:
            if candidate in record_files:
                return record_files[candidate]
        return None

    def download_Database(self, data_dir, symbols, dabase_type, organ=None):
        """Download the reference database via direct Zenodo file links without using the Zenodo API."""
        organ = organ or self.organ
        config = DATABASE_CONFIG[(organ, dabase_type)]
        record_id = config["record_id"]
        base_url = f"https://zenodo.org/records/{record_id}/files"
        def make_url(remote_name):
            return f"{base_url}/{quote(remote_name)}?download=1"
        def download_exact(remote_name, local_name=None):
            return self.download_single_data(
                make_url(remote_name),
                data_dir,
                local_name or remote_name,
            )
        def download_first_available(remote_names, local_name):
            """Support both xxx.topDEGs.csv and xxx_topDEGs.csv naming conventions."""
            errors = []
            for remote_name in remote_names:
                try:
                    return download_exact(remote_name, local_name)
                except RuntimeError as exc:
                    errors.append(str(exc))
            raise RuntimeError(
                f"Unable to download {local_name}; attempted: {', '.join(remote_names)}"
            )
        for symbol in symbols:
            print(f"Downloading reference dataset: {symbol}")
            download_exact(
                f"{symbol}.all.txt",
                f"{symbol}.all.txt",
            )
            download_exact(
                f"{symbol}.Celltype.txt",
                f"{symbol}.Celltype.txt",
            )
            download_first_available(
                [
                    f"{symbol}.topDEGs.csv",
                    f"{symbol}_topDEGs.csv",
                ],
                f"{symbol}.topDEGs.csv",
            )
        orthogroups = download_exact(config["orthogroups"])
        hierarchy = download_exact(config["celltype"])
        return orthogroups, hierarchy

    def map_and_group_by_og(self,df, gene_to_og):
        SYMBOLS = list("ABCDEFGHIJKLMNOPQRSTUVWXYabcdefghijklmnopqrstuvwxy")
        SYMBOL_TO_NUM = {symbol: number for number, symbol in enumerate(SYMBOLS)}
        valid_columns = df.columns[df.columns.isin(gene_to_og)]
        if len(valid_columns) == 0:
            return pd.DataFrame(index=df.index)
        df_num = df.loc[:, valid_columns].apply(lambda column: column.map(SYMBOL_TO_NUM))
        if df_num.isna().any().any():
            invalid = pd.unique(df.loc[:, valid_columns].astype(str).stack())
            invalid = [value for value in invalid if value not in SYMBOL_TO_NUM]
            raise ValueError(
                "The symbolic expression matrix contains values outside A-Y/a-y: " + ", ".join(invalid[:10])
            )
        df_num.columns = [gene_to_og[gene] for gene in valid_columns]
        df_og = (df_num.T.groupby(level=0, sort=True).median().T.round().astype(np.int8))
        return df_og

    def encode_df(self,df, colnames):
        return df.loc[:, colnames].to_numpy(dtype=np.int8,copy=False)

    def normalize(self, s):
        return str(s).strip()

    def fdr_bh(self,pvalues):
        """Benjamini-Hochberg multiple-testing correction."""
        pvalues = np.asarray(pvalues, dtype=float)
        n = len(pvalues)
        order = np.argsort(pvalues)
        ranked = pvalues[order] * n / np.arange(1, n + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        qvalues = np.empty(n)
        qvalues[order] = np.clip(ranked, 0, 1)
        return qvalues

    def evaluate_assignment(self,result_df,score_df,max_shared_og,novelty_cutoff,margin_cutoff=30,alpha=0.05,vote_table=None):
        """Calculate voting significance for each annotation round and classify it as Assigned, De_novo, or Ambiguous."""
        records = []
        if vote_table is None:
            vote_table = pd.crosstab(result_df["query_celltype"],result_df["dataset_celltype"])
        total_all_votes = int(vote_table.to_numpy().sum())
        for query_type in score_df.index:
            scores = score_df.loc[query_type].fillna(0)
            sorted_scores = scores.sort_values(ascending=False)
            if sorted_scores.empty:
                raise ValueError(f"Cell type {query_type} has no reference scores available for annotation")
            assigned_type = sorted_scores.index[0]
            top1_score = sorted_scores.iloc[0]
            top2_score = (sorted_scores.iloc[1]if len(sorted_scores) > 1 else 0)
            margin = top1_score - top2_score
            # Fisher's exact test: votes matching the current query type to the candidate type.
            assigned_votes = int(
                vote_table.at[query_type, assigned_type]
                if query_type in vote_table.index and assigned_type in vote_table.columns
                else 0
            )
            total_votes = int(vote_table.loc[query_type].sum()) if query_type in vote_table.index else 0
            a = assigned_votes
            b = total_votes - assigned_votes
            c = int(vote_table[assigned_type].sum()) - a # Votes from other query types matching this candidate type.
            d = total_all_votes - a - b - c  # Votes from other query types matching other candidate types.
            pvalue = fisher_exact([[a, b], [c, max(0, d)]], alternative="greater").pvalue
            records.append({
                "Query_CellType": query_type,
                "Assigned_CellType": assigned_type,
                "Top1_Score": top1_score,
                "Top2_Score": top2_score,
                "Margin": margin,
                "Assigned_Votes": assigned_votes,
                "Total_Votes": total_votes,
                "Vote_Percent": assigned_votes / total_votes * 100 if total_votes else 0.0,
                "PValue": pvalue,
                "Max_Shared_OG": max_shared_og.get(query_type, 0)
            })
        significance = pd.DataFrame(records)
        significance["QValue"] = self.fdr_bh(significance["PValue"])
        significance["Status"] = np.where(
                significance["Max_Shared_OG"] < novelty_cutoff,
                "De_novo",
                np.where(
                    (significance["QValue"] <= alpha) &
                    (significance["Margin"] >= margin_cutoff),
                    "Assigned",
                    "Ambiguous"
                )
            )
        return significance.set_index("Query_CellType")

    @staticmethod
    @njit(parallel=True, cache=True)
    def top1_l1_parallel(query_encoded, dataset_encoded):
        query_number = query_encoded.shape[0]
        dataset_number = dataset_encoded.shape[0]
        feature_number = query_encoded.shape[1]
        best_distances = np.empty(query_number, dtype=np.int64)
        best_positions = np.empty(query_number, dtype=np.int64)
        for query_position in prange(query_number):
            best_distance = np.int64(9223372036854775807)
            best_position = 0
            for dataset_position in range(dataset_number):
                distance = np.int64(0)
                for feature_position in range(feature_number):
                    difference = (np.int64(dataset_encoded[dataset_position,feature_position])
                                - np.int64(query_encoded[query_position,feature_position]))
                    if difference < 0:
                        difference = -difference
                    distance += difference
                    if distance >= best_distance:
                        break
                if distance < best_distance:
                    best_distance = distance
                    best_position = dataset_position
            best_distances[query_position] = best_distance
            best_positions[query_position] = best_position
        return best_distances, best_positions

    def cal_similarity_fast(self,query, dataset):
        if dataset.empty:
            raise ValueError("Reference dataset is empty")
        common_col = (query.columns.intersection(dataset.columns).sort_values())
        if len(common_col) == 0:
            raise ValueError("The query and dataset have no shared orthogroups")
        # print(f"Number of common OGs: {len(common_col)}")
        query_encoded = np.ascontiguousarray(query.loc[:, common_col].to_numpy(dtype=np.int8,copy=False))
        dataset_encoded = np.ascontiguousarray(dataset.loc[:, common_col].to_numpy(dtype=np.int8,copy=False))
        best_distances, best_positions = self.top1_l1_parallel(query_encoded,dataset_encoded)
        return pd.DataFrame({"query_cell": query.index,
                            "dataset_cell": (dataset.index.to_numpy()[best_positions]),
                            "similarity": (best_distances.astype(np.float64)/ len(common_col))
                            })

    def draw_cell_count(self,data,filename):
        custom_cmap = LinearSegmentedColormap.from_list("custom_blue_gradient", ["#fbf7f7", "#a40303"])
        n_rows, n_cols = data.shape
        cell_size = 0.5 
        figsize = (cell_size * n_cols + 4, cell_size * n_rows + 4)  
        g = sns.clustermap(
            data,
            cmap=custom_cmap,
            annot=True,
            fmt=".1f",
            linewidths=2,
            linecolor="#c9d6df",
            cbar_kws={"shrink": 0.7},
            annot_kws={"size": 16, "color": "black"},
            figsize=figsize,
            dendrogram_ratio=(0.1, 0.1),
            cbar_pos=(1.1, 0.3, 0.03, 0.4), 
            tree_kws={'linewidths': 2, 'colors': '#524748'} ,
            row_cluster=False,
            col_cluster=False,
        )
        g.ax_heatmap.set_title("Top1-predicted Cell Count",fontsize=18,weight="bold",pad=20)
        g.ax_heatmap.set_xlabel("Reference Celltype", fontsize=18,  labelpad=10)
        g.ax_heatmap.set_ylabel("Query Celltype", fontsize=18,  labelpad=10)
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=16)
        plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=16)
        cbar = g.ax_cbar
        cbar.set_ylabel("Cell Count", fontsize=18)
        cbar.tick_params(labelsize=14)
        reordered_data = data
        for y_index, row in enumerate(reordered_data.values):
            max_col_index = row.argmax()
            rect = Rectangle((max_col_index, y_index), 1, 1, fill=False, edgecolor="#fd1158", linewidth=3)
            g.ax_heatmap.add_patch(rect)
        plt.tight_layout()
        plt.savefig(filename, bbox_inches="tight", dpi=300)
        plt.close()

    def draw_cell_percent(self,data,filename):
        custom_cmap = LinearSegmentedColormap.from_list("custom_blue_gradient", ["#ffffff", "#dbd65c", "#5614b0"])
        n_rows, n_cols = data.shape
        cell_size = 0.5 
        figsize = (cell_size * n_cols + 4, cell_size * n_rows + 4)  
        g = sns.clustermap(
            data,
            cmap=custom_cmap,
            annot=True,
            fmt=".1f",
            linewidths=2,
            linecolor="#c9d6df",
            cbar_kws={"shrink": 0.7},
            annot_kws={"size": 14, "color": "black"},
            figsize=figsize,
            dendrogram_ratio=(0.1, 0.1),
            cbar_pos=(1.1, 0.3, 0.03, 0.4), 
            tree_kws={'linewidths': 2, 'colors': '#524748'} ,
            row_cluster=False,
            col_cluster=False,
        )
        g.ax_heatmap.set_title("Prediction of Cell Types (%)",fontsize=18,weight="bold",pad=20)
        g.ax_heatmap.set_xlabel("Reference Celltype", fontsize=18,  labelpad=10)
        g.ax_heatmap.set_ylabel("Query Celltype", fontsize=18,  labelpad=10)
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=16)
        plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=16)
        cbar = g.ax_cbar
        cbar.set_ylabel("Proportion of Top1-predicted Cells (%)", fontsize=16, labelpad=15)
        cbar.tick_params(labelsize=14)
        reordered_data = data
        for y_index, row in enumerate(reordered_data.values):
            max_col_index = row.argmax()
            rect = Rectangle((max_col_index, y_index), 1, 1, fill=False, edgecolor="#fd1158", linewidth=3)
            g.ax_heatmap.add_patch(rect)
        plt.tight_layout()
        plt.savefig(filename, bbox_inches="tight", dpi=300)
        plt.close()

    def draw_cell_percent_clustering(self,data, filename):
        custom_cmap = LinearSegmentedColormap.from_list("custom_blue_gradient", ["#ffffff", "#5b9ff1", "#9c1019"])
        n_rows, n_cols = data.shape
        cell_size = 0.5 
        figsize = (cell_size * n_cols + 4, cell_size * n_rows + 4)  
        g = sns.clustermap(
            data,
            cmap=custom_cmap,
            annot=True,
            fmt=".1f",
            linewidths=2,
            linecolor="#c9d6df",
            cbar_kws={"shrink": 0.7},
            annot_kws={"size": 14, "color": "black"},
            figsize=figsize,
            dendrogram_ratio=(0.1, 0.1),
            cbar_pos=(1.1, 0.3, 0.03, 0.4), 
            tree_kws={'linewidths': 2, 'colors': '#524748'} 
        )
        g.ax_heatmap.set_title("Prediction of Cell Types (%)",fontsize=18,weight="bold",pad=60)
        g.ax_heatmap.set_xlabel("Reference Celltype", fontsize=18,  labelpad=10)
        g.ax_heatmap.set_ylabel("Query Celltype", fontsize=18,  labelpad=10)
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=16)
        plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=16)
        cbar = g.ax_cbar
        cbar.set_ylabel("Proportion of Top1-predicted Cells (%)", fontsize=16, labelpad=15)
        cbar.tick_params(labelsize=14)
        reordered_data = data.iloc[g.dendrogram_row.reordered_ind, g.dendrogram_col.reordered_ind]
        for y_index, row in enumerate(reordered_data.values):
            max_col_index = row.argmax()
            rect = Rectangle((max_col_index, y_index), 1, 1, fill=False, edgecolor="#fd1158", linewidth=3)
            g.ax_heatmap.add_patch(rect)
        plt.tight_layout()
        plt.savefig(filename, bbox_inches="tight", dpi=300)
        plt.close()

    def row_percent(self,data):
        """Convert the matrix to row percentages while retaining all-zero rows as zero."""
        denominator = data.sum(axis=1).replace(0, np.nan)
        return data.div(denominator, axis=0).mul(100).fillna(0)

    def draw_ctype_bar(self,df,filename):
        n_rows, n_cols = df.shape
        custom_colors = ["#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277", "#52AF43", 
                            "#6F9E4C", "#DD9A88", "#F16667", "#E42022", "#F06C45", "#FDBB69", 
                            "#FE982C", "#F78620", "#D9A295", "#B294C7", "#7D54A5", "#9E8099", 
                            "#F0EB99", "#DBB466", "#B15928", "#E58606", "#9A7664", "#5B71AF", 
                            "#559EA7", "#60BE90", "#87C55C", "#A8A965", "#C4709F", "#886A94", 
                            "#2C776F", "#7E8F43", "#D1A323", "#739480", "#3980BE", "#605FAA", 
                            "#93538D", "#D55F67", "#E1556C", "#CF3E88", "#BA6C92", "#A5AA99",
                            "#E58606", "#967568", "#5B76AE", "#54A6A6", "#68C085", "#91C74E",
                            "#B19777", "#C362AC", "#617085", "#4A825A", "#B39B2C", "#9A9B59", 
                            "#388BBB", "#546AB0", "#825097", "#C75D6F", "#E45867", "#D14085", 
                            "#BB6992", "#A5AA99"]
        if n_cols <= 20:
            colors = ["#A6CEE3", "#579CC7", "#3688AD", "#8BC395", "#89CB6C", "#40A635", "#919D5F",
                       "#F99392", "#EB494A", "#E83C2D", "#F79C5D", "#FDA746", "#FE8205", "#E39970", 
                       "#BFA5CF", "#8861AC", "#917099", "#E7E099", "#DEB969", "#B15928"]
        else:
            unique_colors = list(dict.fromkeys(custom_colors))
            if n_cols <= len(unique_colors):
                colors = random.sample(unique_colors, n_cols)
            else:
                colors = sns.color_palette("husl", n_colors=n_cols)
        color_dict = dict(zip(df.columns, colors))
        fig, ax = plt.subplots(figsize=(max(8, n_rows * 0.8), 6))
        bar_width = 0.8
        indices = np.arange(n_rows)
        for i, idx in enumerate(df.index):
            row = df.loc[idx]
            sorted_items = row.sort_values(ascending=True)
            cum_height = 0
            for celltype, val in sorted_items.items():
                ax.bar(i, val, bottom=cum_height, width=bar_width, color=color_dict[celltype], edgecolor='black', linewidth=0.4)
                cum_height += val
        for spine in ax.spines.values():
            spine.set_linewidth(1.8)
        ax.tick_params(axis='both', width=1.6, length=6, labelsize=12)
        ax.set_xticks(indices)
        ax.set_xticklabels(df.index, rotation=45, ha='right', fontsize=12)
        ax.set_ylabel("Proportion of predicted cells(%)", fontsize=14)
        ax.set_xlabel("Query Cell Type", fontsize=14)
        ax.set_title("Assignment Proportion(%)", fontsize=16)
        handles = [plt.Rectangle((0,0),1,1, color=color_dict[ct]) for ct in df.columns]
        ax.legend(handles, df.columns, title="Predicted Cell Type", bbox_to_anchor=(1.02,1), loc='upper left', fontsize=12, title_fontsize=14)
        sns.despine()
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight') 
        plt.close()

    def make_round_output_scores(self,score_df, significance_df, output_types):
        """Generate final summary scores and store de novo/Ambiguous calls in a separate result column."""
        output = score_df.loc[output_types].copy()
        output["de novo"] = 0.0
        output["Ambiguous"] = 0.0
        for celltype in output_types:
            status = significance_df.loc[celltype, "Status"]
            if status in {"De_novo", "Ambiguous"}:
                output.loc[celltype, :] = 0.0
                output.loc[celltype, "de novo" if status == "De_novo" else "Ambiguous"] = 100.0
        return output

    def ancestor_distance(self, child, ancestor, parent_graph=None):
        """Calculate how many hierarchy levels separate an ancestor from a child."""
        parent_graph = self.parent_graph if parent_graph is None else parent_graph
        child = self.normalize(child)
        ancestor = self.normalize(ancestor)
        queue = deque([(child, 0)])
        visited = {child}
        while queue:
            current, distance = queue.popleft()
            if current == ancestor:
                return distance
            for parent in parent_graph.get(current, []):
                if parent not in visited:
                    visited.add(parent)
                    queue.append((parent, distance + 1))
        return None

    def get_prediction_type(self,row):
        status = self.normalize(row["Status"]).lower()
        query_type = self.normalize(row["Query_CellType"])
        assigned_type = self.normalize(row["Assigned_CellType"])
        if status in ["de_novo", "de novo", "denovo"]:
            return "De novo"
        # Exact match.
        if query_type == assigned_type:
            return "Exact"
        # Assigned is an ancestor of Query.
        distance = self.ancestor_distance(query_type, assigned_type)
        if distance is not None:
            return  f"Ancestor_level_{distance}"
        # Assigned is a descendant of Query.
        reverse_distance = self.ancestor_distance(assigned_type, query_type)
        if reverse_distance is not None:
            return "Sublevel"
        return "Cross branch"

    def hierarchy_score(self,prediction_type):
        RHO = 0.95 # rho
        value = str(prediction_type).strip().lower()
        if value in ["exact", "sublevel"]:
            return 1.0
        if value.startswith("ancestor_level_"):
            distance = int(value.split("_")[-1])
            return RHO ** distance
        return 0.0

    def plot_prediction_type(self,final_result, output_file):
        """Plot the C-OHAS score and prediction_type for each query cell type."""
        plot_data = final_result.copy()
        if plot_data.empty:
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.text(0.5, 0.5, "No non-de-novo cell types", ha="center", va="center")
            ax.axis("off")
            fig.tight_layout()
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            plt.close(fig)
            return
        # Convert C-OHAS to numeric; De_novo or "/" becomes NaN.
        plot_data["C_OHAS_Type_Score"] = pd.to_numeric(plot_data["C_OHAS_Type_Score"], errors="coerce").fillna(0)
        # Assign colors according to prediction_type.
        def get_color(prediction_type):
            prediction_type = str(prediction_type)
            if prediction_type == "Exact":
                return "#2a9d8f"
            if prediction_type.startswith("Ancestor_level_"):
                return "#f4a261"
            if prediction_type == "Sublevel":
                return "#457b9d"
            if prediction_type == "De novo":
                return "#8338ec"
            if prediction_type in ["Cross branch", "Unrelated"]:
                return "#e63946"
            return "#999999"
        plot_data = plot_data.sort_values(["C_OHAS_Type_Score", "Query_CellType"],ascending=[True, True])
        colors = plot_data["prediction_type"].apply(get_color)
        fig, ax = plt.subplots(figsize=(7.5, max(5, 0.45 * len(plot_data))))
        ax.barh(plot_data["Query_CellType"],plot_data["C_OHAS_Type_Score"],
                color=colors,edgecolor="black",linewidth=0.5)
        # Add predicted cell types and prediction_type values.
        for position, (_, row) in enumerate(plot_data.iterrows()):
            label = (f'{row["Assigned_CellType"]} | 'f'{row["prediction_type"]}')
            ax.text(min(row["C_OHAS_Type_Score"] + 0.015, 1.02),position,label,va="center",fontsize=9)
        ax.set_xlim(0, 1.35)
        ax.set_xlabel("C-OHAS score of each celltype")
        ax.set_ylabel("Query cell type")
        ax.set_title("C-OHAS Score and Predicted Cell Type")
        ax.axvline(0.8,linestyle="--",color="#666666",linewidth=1)
        sns.despine()
        fig.tight_layout()
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def Genarate_Cell_Seq(self, input_q, symbol_q, out_data, filter_list):
        """Convert a user-provided h5ad file into CellBLASTer symbolic, cell-type, and DEG tables."""
        os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
        try:
            import scanpy as sc
        except ImportError as exc:
            raise ImportError("scanpy is missing; install the CellBLASTer runtime dependencies first") from exc

        input_q = Path(input_q)
        out_data = Path(out_data)
        out_data.mkdir(parents=True, exist_ok=True)
        if not input_q.is_file():
            raise FileNotFoundError(f"Query data does not exist: {input_q}")

        print(f"01. Reading query data: {input_q}")
        all_txt_path = out_data / f"{symbol_q}.all.txt"
        celltype_txt_path = out_data / f"{symbol_q}.Celltype.txt"
        degs_csv_path = out_data / f"{symbol_q}.topDEGs.csv"
        adata = sc.read_h5ad(input_q)
        if "Celltype" not in adata.obs.columns:
            raise ValueError("The input h5ad adata.obs must contain a 'Celltype' column")
        if adata.obs["Celltype"].isna().any():
            raise ValueError("The 'Celltype' column in the input h5ad must not contain missing values")
        adata.obs["Celltype"] = adata.obs["Celltype"].astype(str)
        if adata.obs_names.has_duplicates or adata.var_names.has_duplicates:
            raise ValueError("Cell names and gene names in the input h5ad must be unique")

        adata.var["name"] = adata.var_names.astype(str)
        if filter_list:
            filter_pattern = "|".join(re.escape(item) for item in filter_list if item)
            if filter_pattern:
                print(f"Filtering genes containing the following keywords: {', '.join(filter_list)}")
                filtered_genes = adata.var["name"].str.contains(
                    filter_pattern, case=False, na=False, regex=True
                )
                adata = adata[:, ~filtered_genes].copy()

        sc.pp.filter_genes(adata, min_cells=5)
        if adata.n_vars == 0:
            raise ValueError("No genes remain after filtering")
        if adata.n_obs < 3 or adata.n_vars < 2:
            raise ValueError("The query requires at least 3 cells and 2 genes that pass filtering")
        adata.uns['raw'] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)

        print("02. Running DBSCAN + LOF outlier-cell detection...")
        sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars))
        adata.raw = adata.copy()
        adata = adata[:, adata.var["highly_variable"]]
        sc.pp.scale(adata)
        sc.tl.pca(adata, svd_solver="arpack")

        eps = 10 # First-pass DBSCAN filtering.
        min_samples = 10
        outliers = []
        lof_n_neighbors = 20
        lof_contamination = 0.05
        for cell_type in adata.obs['Celltype'].unique():
            ct_mask = adata.obs['Celltype'] == cell_type
            ct_data = adata[ct_mask]
            X_scaled = StandardScaler().fit_transform(ct_data.obsm['X_pca'])
            db_labels  = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(X_scaled)
            db_outliers = db_labels  == -1
            lof_labels = LocalOutlierFactor(n_neighbors=lof_n_neighbors, contamination=lof_contamination).fit_predict(X_scaled)
            lof_outliers = lof_labels == -1
            combined_outliers = db_outliers | lof_outliers
            ct_indices = np.where(ct_mask)[0]
            outliers.extend(ct_indices[combined_outliers])

        outlier_mask = np.zeros(adata.n_obs, dtype=bool)
        outlier_mask[outliers] = True
        adata.obs['outlier'] = 'Inlier'
        adata.obs.loc[outlier_mask, 'outlier'] = 'Outlier'
        # Update the adata object.
        adata = adata.raw.to_adata()
        adata.X = adata.uns['raw'].copy()
        adata = adata[ adata.obs['outlier'] == 'Inlier'].copy()
        adata.raw = adata.copy()
        sc.pp.normalize_total(adata,target_sum=10000)
        sc.pp.log1p(adata)

        print("03. Running Wilcoxon DEG analysis...")
        sc.tl.rank_genes_groups(adata, groupby='Celltype',  method='wilcoxon', key='rank_genes_groups')
        degs = sc.get.rank_genes_groups_df(adata, group = None, pval_cutoff=0.05, log2fc_min=0.5) 
        top200_degs = degs.groupby('group').head(150) 
        top200_degs.to_csv(degs_csv_path, sep="\t")

        gene_3 = top200_degs['names'].unique().tolist()
        # Replace with the new adata object.
        adata = adata.raw.to_adata()
        adata1 = adata[:, gene_3].copy()
        sc.pp.normalize_total(adata1,target_sum=10000)
        sc.pp.log1p(adata1)
        sc.pp.scale(adata1)

        print("04. Generating the symbolic expression matrix...")
        df = pd.DataFrame(adata1.X,index = adata1.obs.index,columns = adata1.var.index)
        df = df.rank(pct=True, axis=1) * 100
        df = df.round().astype(int)
        symbols = [chr(ord('A') + i) for i in range(25)] + [chr(ord('a') + i) for i in range(25)]
        bins = np.linspace(0, 100, 51) 
        df_symbolic = df.apply(lambda col: pd.cut(col, bins=bins, labels=symbols, include_lowest=True))
        # Save as the wheat percentile file.
        df_symbolic.to_csv(all_txt_path, sep="\t")
        Celltype = adata1.obs[['Celltype']]
        Celltype['Cell'] = df.index
        Celltype.to_csv(celltype_txt_path, sep="\t", index=False)
        return {
            "all": str(all_txt_path),
            "celltype": str(celltype_txt_path),
            "degs": str(degs_csv_path),
        }

    def load_databse(self, data_dir, gene_to_og, symbols):
        """Read and harmonize reference datasets; retain the legacy method spelling for compatibility."""
        data_dir = Path(data_dir)
        datasets = []
        Celltype_data = []
        all_DEG = []
        for i, symbol in enumerate(symbols):
            prefix = f"d{i+1}"
            all_file = data_dir / f"{symbol}.all.txt"
            celltype_file = data_dir / f"{symbol}.Celltype.txt"
            deg_file = data_dir / f"{symbol}.topDEGs.csv"
            missing = [str(path) for path in (all_file, celltype_file, deg_file) if not path.is_file()]
            if missing:
                raise FileNotFoundError("Missing reference data files: " + ", ".join(missing))
            print(f"Processing {symbol} with prefix {prefix}...")
            ds = pd.read_csv(all_file, sep='\t', index_col=0)
            ct = pd.read_csv(celltype_file, sep='\t')
            deg = pd.read_csv(deg_file, sep='\t',index_col=0)
            if not {"Cell", "Celltype"}.issubset(ct.columns):
                raise ValueError(f"{celltype_file} must contain the Cell and Celltype columns")
            if not {"group", "names"}.issubset(deg.columns):
                raise ValueError(f"{deg_file} must contain the group and names columns")
            ds.columns = ds.columns.astype(str)
            ct["Cell"] = ct["Cell"].astype(str)
            ct["Celltype"] = ct["Celltype"].astype(str)
            deg["group"] = deg["group"].astype(str)
            ds.index = [f"{prefix}_{str(idx)}" for idx in ds.index]
            ds = self.map_and_group_by_og(ds, gene_to_og)
            if ds.empty or ds.shape[1] == 0:
                raise ValueError(f"{symbol} has no genes that map to the selected orthogroups")
            ct['Cell'] = [f"{prefix}_{str(idx)}" for idx in ct['Cell']]
            datasets.append(ds)
            Celltype_data.append(ct)
            all_DEG.append(deg)
        if not datasets:
            raise ValueError("No reference datasets were loaded successfully")
        Celltype_data = pd.concat(Celltype_data, ignore_index=True)
        print(f"Merging complete! Total cells: {sum(len(dataset) for dataset in datasets)}")
        return Celltype_data,datasets,all_DEG

    def load_database(self, data_dir, gene_to_og, symbols):
        return self.load_databse(data_dir, gene_to_og, symbols)

    def Generate_User_Reference(
        self,
        reference_adata,
        reference_symbol,
        data_dir,
        filter_list,
    ):
        """
        Convert a user-provided reference h5ad file into the three
        CellBLASTer reference database files.
        """
        if reference_adata is None and reference_symbol is None:
            return None
        if reference_adata is None or reference_symbol is None:
            raise ValueError(
                "reference_adata and reference_symbol must be provided together"
            )
        reference_path = Path(reference_adata).expanduser().resolve()
        if not reference_path.is_file():
            raise FileNotFoundError(f"reference h5ad file does not exist：{reference_path}")
        if not re.fullmatch(
            r"[A-Za-z0-9][A-Za-z0-9._-]*",
            reference_symbol,
        ):
            raise ValueError("reference symbol can only contain letters, numbers, dots, underscores and hyphens")
        data_dir = Path(data_dir)
        data_dir.mkdir(parents=True, exist_ok=True)
        print(f"Generate a user-defined reference dataset："f"{reference_symbol}")
        self.Genarate_Cell_Seq(
            input_q=reference_path,
            symbol_q=reference_symbol,
            out_data=data_dir,
            filter_list=filter_list,
        )
        print(f"User-defined reference dataset has been written：{data_dir}")
        return reference_symbol


    def Annotation(
        self,
        output_path=None,
        symbols=None,
        dabase_type=None,
        query=None,
        query_symbol=None,
        filter_keywords=None,
        organ=None,
        sample_ratio=None,
    ):
        """Run downloads, query preprocessing, iterative annotation, and result visualization."""
        output_path = Path(output_path or self.output_path).expanduser().resolve()
        symbols = list(symbols or self.symbols)
        database_type = dabase_type or self.database_type
        query_path = Path(query or self.query).expanduser().resolve()
        query_symbol = query_symbol or self.query_symbol
        filter_keywords = self.filter_keywords if filter_keywords is None else filter_keywords
        organ = organ or self.organ
        reference_adata = self.reference_adata
        reference_symbol = self.reference_symbol

        if not symbols:
            raise ValueError("At least one reference dataset symbol must be specified")
        if (organ, database_type) not in DATABASE_CONFIG:
            raise ValueError("organ must be Root/Leaf/Flower, and database_type must be Dicot/Monocot")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", query_symbol):
            raise ValueError("query_symbol may contain only letters, numbers, periods, underscores, and hyphens")
        if not query_path.is_file():
            raise FileNotFoundError(f"Query h5ad file does not exist: {query_path}")
        output_path.mkdir(parents=True, exist_ok=True)

        data_dir = output_path / "01.DataBase"
        data_dir.mkdir(parents=True, exist_ok=True)
        print(f"Downloading from the {organ}/{database_type} database: {len(symbols)} reference dataset(s) to {data_dir}...")
        orthogroups_path, hierarchy_path = self.download_Database(data_dir, symbols, database_type, organ=organ)
        all_reference_symbols = list(symbols)
        user_reference_symbol = self.Generate_User_Reference(
            reference_adata=reference_adata,
            reference_symbol=reference_symbol,
            data_dir=data_dir,
            filter_list=filter_keywords,
        )
        if user_reference_symbol is not None:all_reference_symbols.append(user_reference_symbol)
        print("All reference datasets involved in the annotation:",all_reference_symbols,)

        print("Reading orthogroup information...")
        records = []
        with open(orthogroups_path, "r", encoding="utf-8") as handle:
            for line in handle:
                if ":" not in line:
                    continue
                og, genes_text = line.strip().split(":", 1)
                records.extend((og, gene) for gene in genes_text.strip().split())
        if not records:
            raise ValueError(f"The orthogroups file is empty or cannot be parsed: {orthogroups_path}")
        og_table = pd.DataFrame(records, columns=["OG", "Gene"])
        og_table["Gene"] = og_table["Gene"].str.replace("gene=", "", regex=False)
        og_table = og_table.drop_duplicates(subset="Gene", keep="first")
        gene_to_og = dict(zip(og_table["Gene"], og_table["OG"]))

        print("Loading the reference database...")
        celltype_data, datasets, all_deg_tables = self.load_database(data_dir, gene_to_og, all_reference_symbols,)
        query_dir = output_path / "02.QueryData"
        print(f"Generating query intermediate files in {query_dir}...")
        self.Genarate_Cell_Seq(
            query_path, query_symbol, query_dir, filter_keywords
        )

        # Store all iterative results under 03.Blast_Result/<query_symbol>_blast.
        # round_0 stores shared OGs and preliminary de novo calls; round_1... stores each iteration.
        result_root = output_path / "03.Blast_Result"
        result_root.mkdir(parents=True, exist_ok=True)
        round0_dir = result_root / "round_0"
        round0_dir.mkdir(parents=True, exist_ok=True)
        print(f"Starting annotation; results will be written to {result_root}...")

        all_deg = pd.concat(all_deg_tables, ignore_index=True)
        # all_deg["group"] = all_deg["group"].astype(str)
        all_deg["names"] = all_deg["names"].map(gene_to_og)
        query_deg = pd.read_csv(
            query_dir / f"{query_symbol}.topDEGs.csv", sep="\t", index_col=0
        )
        if not {"group", "names"}.issubset(query_deg.columns):
            raise ValueError("The query DEG file must contain the group and names columns")
        # query_deg["group"] = query_deg["group"].astype(str)
        query_deg["names"] = query_deg["names"].map(gene_to_og)

        query_all = pd.read_csv(
            query_dir / f"{query_symbol}.all.txt", sep="\t", index_col=0
        )
        celltype_query = pd.read_csv(
            query_dir / f"{query_symbol}.Celltype.txt", sep="\t"
        )
        if not {"Cell", "Celltype"}.issubset(celltype_query.columns):
            raise ValueError("The query Celltype file must contain the Cell and Celltype columns")
        query_all.index = query_all.index.astype(str)
        query_all.columns = query_all.columns.astype(str)
        celltype_query["Cell"] = celltype_query["Cell"].astype(str)
        celltype_query["Celltype"] = celltype_query["Celltype"].astype(str)
        query_all = self.map_and_group_by_og(query_all, gene_to_og)
        if query_all.empty or query_all.shape[1] == 0:
            raise ValueError("The query has no genes that map to the selected orthogroups")

        # Retain this variable for subsequent iterative annotation.
        all_query_types = list(pd.unique(celltype_query["Celltype"].dropna()))
        # Shared-OG calculation identical to blast_sig_new.py.
        shared_og_all = pd.DataFrame(0,index=query_deg["group"].unique(),columns=all_deg["group"].unique(),)
        for g1 in query_deg["group"].unique():
            query_ogs = {
                x
                for x in query_deg.loc[
                    query_deg["group"] == g1,
                    "names",
                ]
                if pd.notna(x)
            }
            for g2 in all_deg["group"].unique():
                reference_ogs = {
                    x
                    for x in all_deg.loc[
                        all_deg["group"] == g2,
                        "names",
                    ]
                    if pd.notna(x)
                }
                shared = query_ogs & reference_ogs
                shared_og_all.loc[g1, g2] = len(shared)
        shared_og_all = (shared_og_all.sort_index().sort_index(axis=1))
        shared_og_all.to_csv(round0_dir / "Shared_OMG_all_CType.csv")

        datasets = [dataset.loc[:, sorted(query_all.columns.intersection(dataset.columns))] for dataset in datasets]
        datasets = [dataset for dataset in datasets if dataset.shape[1] > 0]
        if not datasets:
            raise ValueError("The query shares no orthogroups with any reference dataset")

        if sample_ratio is not None:
            if not 0 < sample_ratio < 1:
                raise ValueError("sample_ratio must be a value between 0 and 1")
            celltype_data = celltype_data.copy()
            celltype_data["Dataset"] = celltype_data["Cell"].str.extract(r"^(d\d+)")
            strata = celltype_data["Dataset"] + "_" + celltype_data["Celltype"].astype(str)
            try:
                celltype_data, _ = train_test_split(
                    celltype_data,
                    train_size=sample_ratio,
                    stratify=strata,
                    random_state=42,
                )
            except ValueError as exc:
                raise ValueError(f"Unable to stratify and downsample reference data by Dataset + Celltype: {exc}") from exc
            sampled_cell_ids = set(celltype_data["Cell"])
            datasets = [
                dataset.loc[dataset.index.isin(sampled_cell_ids)] for dataset in datasets
            ]
            datasets = [dataset for dataset in datasets if not dataset.empty]
            print(f"Reference data stratified and sampled by dataset and cell type: {sample_ratio:.1%}")
        else:
            print("No downsampling; retaining all reference cells")

        reference_datasets = datasets
        reference_celltypes = celltype_data.copy()
        query_celltype_map = dict(
            zip(celltype_query["Cell"], celltype_query["Celltype"])
        )
        reference_celltype_map = dict(
            zip(reference_celltypes["Cell"], reference_celltypes["Celltype"])
        )
        print("Number of query cells:", len(query_all))
        print("Number of query cell clusters:", len(all_query_types))
        print("Reference dataset cell counts:", [len(dataset) for dataset in reference_datasets])

        margin_cutoff = 30
        significance_alpha = 0.05
        resolved_types = []
        remaining_types = all_query_types.copy()
        significance_history = []
        round_score_tables = []

        max_shared_og = shared_og_all.max(axis=1)
        q1 = max_shared_og.quantile(0.25)
        q3 = max_shared_og.quantile(0.75)
        novelty_cutoff = min(15,max(1.0, float(q1 - 1.5 * (q3 - q1))))
        print("Shared-OG threshold for de novo detection:", novelty_cutoff)
        pre_denovo_types = max_shared_og[
            (max_shared_og < novelty_cutoff) & (max_shared_og < 10)
        ].index.tolist()

        if pre_denovo_types:
            print("De novo cell types identified in round 0:", pre_denovo_types)
            pre_denovo_df = pd.DataFrame(
                {
                    "Assigned_CellType": "de novo",
                    "Top1_Score": np.nan,
                    "Top2_Score": np.nan,
                    "Margin": np.nan,
                    "Assigned_Votes": 0,
                    "Total_Votes": 0,
                    "Vote_Percent": 0.0,
                    "PValue": np.nan,
                    "QValue": np.nan,
                    "Max_Shared_OG": max_shared_og.loc[pre_denovo_types],
                    "Status": "De_novo",
                    "Round": 0,
                },
                index=pre_denovo_types,
            )
            pre_denovo_df.index.name = "Query_CellType"
            significance_history.append(pre_denovo_df)
            pre_denovo_df.to_csv(round0_dir / "Assignment_significance_round_0.csv")
            pre_denovo_scores = pd.DataFrame(
                {"de novo": 100.0}, index=pre_denovo_types
            )
            pre_denovo_scores.to_csv(round0_dir / "De_novo_similar_round_0.csv")
            round_score_tables.append(pre_denovo_scores)
            remaining_types = [
                celltype for celltype in remaining_types if celltype not in pre_denovo_types
            ]
        else:
            print("No de novo cell types identified in round 0")

        round_number = 1
        while remaining_types:
            phase = f"round_{round_number}"
            round_dir = result_root / phase
            round_dir.mkdir(parents=True, exist_ok=True)
            print(f"\n========== Round {round_number} ==========")
            print("Cell clusters pending analysis:", remaining_types)

            active_query_meta = celltype_query[
                celltype_query["Celltype"].isin(remaining_types)
            ]
            active_query = query_all.loc[
                query_all.index.isin(active_query_meta["Cell"])
            ]
            if active_query.empty:
                raise ValueError("No query cells are available for analysis in the current round; check the Cell column and matrix index")

            active_reference_meta = reference_celltypes[
                ~reference_celltypes["Celltype"].isin(resolved_types)
            ]
            active_reference_ids = set(active_reference_meta["Cell"])
            active_datasets = [
                dataset.loc[dataset.index.isin(active_reference_ids)]
                for dataset in reference_datasets
            ]
            active_datasets = [dataset for dataset in active_datasets if not dataset.empty]
            if not active_datasets:
                # If all reference types with the same name are excluded, fall back to the full reference set so the workflow can complete.
                active_datasets = reference_datasets

            result_df = pd.concat(
                [
                    self.cal_similarity_fast(active_query, dataset)
                    for dataset in active_datasets
                ],
                ignore_index=True,
            )
            result_df["query_celltype"] = result_df["query_cell"].map(
                query_celltype_map
            )
            result_df["dataset_celltype"] = result_df["dataset_cell"].map(
                reference_celltype_map
            )
            if result_df[["query_celltype", "dataset_celltype"]].isna().any().any():
                raise ValueError("Cell-type mapping failed; check the Cell identifiers in the Celltype file")
            result_df.to_csv(round_dir / f"01.All_cell_match_{phase}.csv", index=False)

            pivot_df = result_df.pivot_table(
                index="query_celltype",
                columns="dataset_celltype",
                aggfunc="size",
                fill_value=0,
            ).astype(int)
            if pivot_df.empty:
                raise ValueError("No valid cell-type matches were generated in the current round")
            pivot_normalized = self.row_percent(pivot_df).round(2)
            pivot_df.to_csv(round_dir / f"02.Top1-predicted_count_{phase}.csv")
            pivot_normalized.to_csv(round_dir / f"03.Top1-predicted_percent_{phase}.csv")
            self.draw_cell_count(pivot_df, round_dir / f"02.Top1-predicted_count_{phase}.png")
            self.draw_cell_percent(
                pivot_normalized, round_dir / f"03.Top1-predicted_percent_{phase}.png"
            )

            shared_og_current = shared_og_all.reindex(
                index=pivot_df.index, columns=pivot_df.columns, fill_value=0
            ).astype(int)
            shared_og_current.to_csv(round_dir / f"04.Shared_OG_{phase}.csv")
            self.draw_cell_count(
                shared_og_current, round_dir / f"04.Shared_OG_{phase}.png"
            )
            weighted_raw = shared_og_current * pivot_normalized
            score_df = self.row_percent(weighted_raw)
            score_df.to_csv(round_dir / f"05.Blast_similarity_{phase}.csv")
            self.draw_cell_percent(
                score_df, round_dir / f"05.Blast_similarity_{phase}.png"
            )

            significance_df = self.evaluate_assignment(
                result_df=result_df,
                score_df=score_df,
                max_shared_og=max_shared_og,
                novelty_cutoff=novelty_cutoff,
                margin_cutoff=margin_cutoff,
                alpha=significance_alpha,
                vote_table=pivot_df,
            )
            significance_df["Round"] = round_number
            assigned_new = significance_df.index[
                significance_df["Status"] == "Assigned"
            ].tolist()
            denovo_new = significance_df.index[
                significance_df["Status"] == "De_novo"
            ].tolist()
            resolved_new = assigned_new + denovo_new
            if not resolved_new:
                ambiguous_candidates = significance_df.index[
                    significance_df["Status"] == "Ambiguous"
                ].tolist()
                if not ambiguous_candidates:
                    raise RuntimeError("No cell types can be advanced in the current round")
                forced_type = significance_df.loc[
                    ambiguous_candidates, "Margin"
                ].idxmax()
                significance_df.loc[forced_type, "Status"] = "Assigned"
                resolved_new = [forced_type]

            significance_df.to_csv(
                round_dir / f"Assignment_significance_{phase}.csv"
            )
            significance_history.append(significance_df.copy())
            print("Cell clusters assigned in this round:", resolved_new)
            round_scores = self.make_round_output_scores(
                score_df, significance_df, resolved_new
            )
            round_scores.to_csv(round_dir / f"06.Assigned_celltype_in_{phase}.csv")
            self.draw_cell_percent(
                round_scores, round_dir / f"06.Assigned_celltype_in_{phase}.png"
            )
            round_score_tables.append(round_scores)
            resolved_types.extend(
                celltype for celltype in resolved_new if celltype not in resolved_types
            )
            remaining_types = [
                celltype for celltype in remaining_types if celltype not in resolved_new
            ]
            round_number += 1

        final_dir = result_root / "final_visual"
        final_dir.mkdir(parents=True, exist_ok=True)
        if not round_score_tables or not significance_history:
            raise RuntimeError("No annotation results were generated for summarization")

        combined_df = pd.concat(round_score_tables, axis=0).fillna(0)
        combined_df = combined_df.groupby(level=0, sort=False).mean()
        combined_df.to_csv(final_dir / "Proportion_Top1-predicted.csv")
        self.draw_cell_percent(combined_df, final_dir / "2.Proportion_Top1-predicted.png")
        self.draw_ctype_bar(combined_df, final_dir / "1.Assigned_Ctype_bar.png")
        self.draw_cell_percent(combined_df, final_dir / "2.Proportion_Top1-predicted.pdf")
        self.draw_ctype_bar(combined_df, final_dir / "1.Assigned_Ctype_bar.pdf")
        if combined_df.shape[0] > 1 and combined_df.shape[1] > 1:
            self.draw_cell_percent_clustering(combined_df, final_dir / "3.Proportion_Top1-predicted_clustered.png")
            self.draw_cell_percent_clustering(combined_df, final_dir / "3.Proportion_Top1-predicted_clustered.pdf")

        final_result = pd.concat(significance_history).sort_values("Round")
        final_result = final_result[
            ~final_result.index.duplicated(keep="last")
        ].copy()
        final_result.index.name = "Query_CellType"
        final_columns = [
            "Query_CellType",
            "Assigned_CellType",
            "Status",
            "Margin",
            "PValue",
            "QValue",
            "Max_Shared_OG",
            "Round",
            "Top1_Score",
            "Top2_Score",
            "Assigned_Votes",
            "Total_Votes",
            "Vote_Percent",
        ]
        final_result = final_result.reset_index().reindex(columns=final_columns)

        celltype_hierarchy = pd.read_csv(hierarchy_path, sep="\t")
        if not {"ct1", "ct2"}.issubset(celltype_hierarchy.columns):
            raise ValueError(f"The hierarchy file must contain the ct1 and ct2 columns: {hierarchy_path}")
        ct_pairs = {
            frozenset(self.normalize(item) for item in pair)
            for pair in celltype_hierarchy[["ct1", "ct2"]].values
        }
        self.parent_graph = defaultdict(set)
        for parent, child in celltype_hierarchy[["ct1", "ct2"]].values:
            self.parent_graph[self.normalize(child)].add(self.normalize(parent))

        denovo_mask = final_result["Status"].eq("De_novo")
        valid_mask = ~denovo_mask
        final_result["Matched"] = pd.Series(
            "/", index=final_result.index, dtype=object
        )
        final_result.loc[valid_mask, "Matched"] = final_result.loc[valid_mask].apply(
            lambda row: self.normalize(row["Query_CellType"])
            == self.normalize(row["Assigned_CellType"])
            or frozenset(
                [
                    self.normalize(row["Query_CellType"]),
                    self.normalize(row["Assigned_CellType"]),
                ]
            )
            in ct_pairs,
            axis=1,
        )
        accuracy = (
            final_result.loc[valid_mask, "Matched"].astype(bool).mean()
            if valid_mask.any()
            else np.nan
        )
        final_result["Accuracy"] = pd.Series(
            "/", index=final_result.index, dtype=object
        )
        final_result.loc[valid_mask, "Accuracy"] = accuracy

        slash_columns = [
            "Assigned_CellType",
            "Margin",
            "Round",
            "Top1_Score",
            "Top2_Score",
            "Assigned_Votes",
            "Total_Votes",
            "Vote_Percent",
        ]
        final_result[slash_columns] = final_result[slash_columns].astype(object)
        final_result.loc[denovo_mask, slash_columns] = "/"
        final_result["prediction_type"] = final_result.apply(
            self.get_prediction_type, axis=1
        )
        final_result["QValue_Significant"] = np.where(
            pd.to_numeric(final_result["QValue"], errors="coerce") < 0.05,
            "Yes",
            "No",
        )
        #annotation_path = final_dir / "Final_annotation.csv"
        #final_result.to_csv(annotation_path, index=False)
        print("Accuracy after excluding de novo calls:", accuracy)

        evaluation_df = final_result[final_result["Status"] != "De_novo"].copy()
        evaluation_df["Hierarchy_Score"] = evaluation_df["prediction_type"].apply(
            self.hierarchy_score
        )
        qvalue = pd.to_numeric(evaluation_df["QValue"], errors="coerce")
        evaluation_df["QValue_Evidence"] = (
            -np.log10(qvalue.clip(lower=1e-300)) / -np.log10(0.05)
        ).clip(lower=0, upper=1)
        margin = pd.to_numeric(evaluation_df["Margin"], errors="coerce")
        evaluation_df["Margin_Evidence"] = margin.div(margin_cutoff).clip(
            lower=0, upper=1
        )
        vote_percent = pd.to_numeric(
            evaluation_df["Vote_Percent"], errors="coerce"
        )
        evaluation_df["Vote_Evidence"] = np.sqrt(
            vote_percent.div(25).clip(lower=0, upper=1)
        )
        evaluation_df["Combined_Evidence"] = (
            0.25 * evaluation_df["QValue_Evidence"]
            + 0.375 * evaluation_df["Vote_Evidence"]
            + 0.375 * evaluation_df["Margin_Evidence"]
        )
        evaluation_df["C_OHAS_Type_Score"] = evaluation_df["Hierarchy_Score"] * (
            0.8 + 0.2 * evaluation_df["Combined_Evidence"]
        )
        c_ohas_closed = evaluation_df["C_OHAS_Type_Score"].mean()
        evaluation_df["Final_C_OHAS"] = c_ohas_closed
        evaluation_path = final_dir / "Final_evaluation.csv"
        evaluation_df.to_csv(evaluation_path, index=False)
        prediction_plot = final_dir / "4.C_OHAS_prediction_type.png"
        self.plot_prediction_type(evaluation_df, prediction_plot)
        self.plot_prediction_type(
            evaluation_df, final_dir / "4.C_OHAS_prediction_type.pdf"
        )
        print("C_OHAS_Type_Score:", c_ohas_closed)
        return {
            "database_dir": str(data_dir),
            "query_dir": str(query_dir),
            "result_dir": str(result_root),
            #"annotation": str(annotation_path),
            "evaluation": str(evaluation_path),
        }

def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Download a selected CellBLASTer reference database, annotate cell types "
            "from an h5ad query, and generate result visualizations."
        )
    )
    parser.add_argument(
        "-t",
        "--database_type",
        "--dabase_type",
        dest="database_type",
        required=True,
        choices=["Dicot", "Monocot"],
        help="Reference database clade: Dicot or Monocot.",
    )
    parser.add_argument(
        "-p",
        "--organ",
        required=True,
        choices=["Root", "Leaf", "Flower"],
        help="Reference organ database: Root, Leaf, or Flower.",
    )
    parser.add_argument(
        "-s",
        "--symbols",
        nargs="+",
        required=True,
        help="Reference dataset symbols, e.g. SRP406470 SRP390780.",
    )
    parser.add_argument(
        "-o", "--output_path", default="./CellBLASTer_output", help="Output directory."
    )
    parser.add_argument(
        "-q", "--query", required=True, help="Path to the input .h5ad file."
    )
    parser.add_argument(
        "-qs", "--query_symbol", required=True, help="Unique query name/prefix."
    )
    parser.add_argument(
        "-f",
        "--filter_keywords",
        nargs="*",
        default=DEFAULT_NONCODING_RNA_KEYWORDS.copy(),
        help=(
            "Case-insensitive non-coding RNA gene-name substrings to remove. "
            "Use -f with no values to disable filtering, or provide custom "
            "keywords to replace the default list."
        ),
    )
    parser.add_argument(
        "--sample-ratio",
        type=float,
        default=None,
        help="Optional stratified reference downsampling ratio between 0 and 1.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=None,
        help="CPU cores used by the Numba similarity calculation (default: min(30, available)).",
    )
    parser.add_argument(
    "--reference-adata",
    default=None,
    help=(
        "Optional path to a user-defined reference h5ad file. "
        "Its adata.obs must contain the Celltype column."
        ),
    )
    parser.add_argument(
        "--reference-symbol",
        default=None,
        help=(
            "Sample name used to generate files for the "
            "user-defined reference dataset."
        ),
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    cellblaster = CellBlaster(
        database_type=args.database_type,
        organ=args.organ,
        symbols=args.symbols,
        output_path=args.output_path,
        query=args.query,
        query_symbol=args.query_symbol,
        filter_keywords=args.filter_keywords,
        reference_adata=args.reference_adata,
        reference_symbol=args.reference_symbol,
        n_jobs=args.n_jobs,
    )
    cellblaster.Annotation(sample_ratio=args.sample_ratio)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
