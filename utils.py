import pandas as pd
import os
import glob
import pypgx
import re
from itertools import product
import plotly.graph_objects as go


def generate_pharmvar_historical_data(pharmvar_dir):
    first_last_pharmvar_version = ["1.1.9.2","6.1.3"]
    df_star_alleles = pd.DataFrame(columns = ["Allele", "Gene", "Version", "Date"])
    df_list = []
    df_list2 = []
    for version_data_dir in os.listdir(pharmvar_dir):
        version_date_path = os.path.join(pharmvar_dir, version_data_dir)
        date = version_data_dir.split('_')[1]

        
        # Check if the directory follows the pharmvar-version_date structure
        if os.path.isdir(version_date_path) and version_data_dir.startswith('pharmvar-'):
            # Iterate over version directories within the version directory
            for pharmvar_version_path in os.listdir(version_date_path):
                pharmvar_path = os.path.join(version_date_path, pharmvar_version_path)
                date = version_data_dir.split("_")[1]
                pharmvar_version = version_data_dir.split("_")[0].split("-")[1]

                if (len(pharmvar_version) == 3 and pharmvar_version[-1] == "0") or pharmvar_version in first_last_pharmvar_version:

                    for gene in os.listdir(pharmvar_path):
                        if str(gene) != "DPYD":
                            gene_dir = os.path.join(pharmvar_path, gene)
                        
                            # Check if the gene directory contains the GRCh38 subdirectory
                            grch38_dir = os.path.join(gene_dir, 'GRCh38')
                            
                            if os.path.isdir(grch38_dir):
                                # I dont know what the file is called, but it's the only tsv file in this directory
                                tsv_file = glob.glob(os.path.join(grch38_dir, '*.tsv'))[0]
                                df = pd.read_csv(tsv_file, sep='\t', header=1)
                                vcf_files = glob.glob(os.path.join(grch38_dir, '*.vcf'))
                                total_number_of_star_alleles = len(vcf_files) + 1 # +1 for the reference allele
                                star_allele_list = df["Haplotype Name"].to_list()
                                major_alleles_raw = [star_allele.split(".")[0] for star_allele in star_allele_list]
                                major_alleles = [major_allele.rstrip(major_allele[-1]) if major_allele[-1].isalpha() else major_allele for major_allele in major_alleles_raw]
                                major_alleles = set(major_alleles)

                                number_of_major_alleles = len(major_alleles)

                                gene = df["Gene"].unique()[0]
                                df_star_alleles["total_star_allele_count"] = total_number_of_star_alleles
                                df_star_alleles["major_allele_count"] = number_of_major_alleles
                                df_star_alleles["Allele"] = df["Haplotype Name"]
                                df_star_alleles["Gene"] = gene  
                                df_star_alleles["Version"] = version_data_dir.split("_")[0].split("-")[1]
                                df_star_alleles["Date"] = date

                                df_allele_version = pd.DataFrame(df_star_alleles, columns = ["Allele", "Version"])
                                df["Version"] = version_data_dir.split("_")[0].split("-")[1]
                                df["Date"] = date
                                df_list2.append(df)

                                # Drop rows with NaN values, drop duplicates
                                df2 = df_star_alleles.dropna()
                                df3 = df2.drop_duplicates(subset = "Allele")
                                ## df.append has been removed, use pd.concat instead: put the dataframes in a list and then concatenate them.
                                df_list.append(df3)
    final_allele_version_df = pd.concat(df_list2, ignore_index = True)
    final_allele_version_df = final_allele_version_df.rename(columns = {"Haplotype Name": "Allele"})
    final_allele_version_df["Date"] = pd.to_datetime(final_allele_version_df["Date"], format="%d%m%Y")
    final_allele_version_df["Date"] = final_allele_version_df["Date"].dt.strftime("%Y-%m-%d")
    final_allele_version_df = final_allele_version_df[["Allele", "Gene", "Date", "Version"]]
    return final_allele_version_df

def keep_major_allele(allele):
    allele = allele.split(".")[0]
    return allele.rstrip(allele[-1]) if allele[-1].isalpha() else allele

def process_diplotype_patterns(diplotype: str) -> str:
    def remove_letters(string):
        return ''.join(char for char in string if not char.isalpha())
    patterns = [
        # a/b(c)
        (r'([A-Za-z0-9]+)/([A-Za-z0-9]+)\(([A-Za-z0-9]+)\)',
         lambda m: f"{remove_letters(m.group(1))}/({remove_letters(m.group(2))},{remove_letters(m.group(3))})"),
        # a/(b,c)
        (r'([A-Za-z0-9]+)/\(([A-Za-z0-9]+),([A-Za-z0-9]+)\)',
         lambda m: f"{remove_letters(m.group(1))}/({remove_letters(m.group(2))},{remove_letters(m.group(3))})"),
               
        # a(b)/c
        (r'([A-Za-z0-9]+)\(([A-Za-z0-9]+)\)/([A-Za-z0-9]+)',
         lambda m: f"({remove_letters(m.group(1))},{remove_letters(m.group(2))})/({remove_letters(m.group(3))})"),
        
        # a/(b)
        (r'([A-Za-z0-9]+)/\(([A-Za-z0-9]+)\)',
         lambda m: f"{remove_letters(m.group(1))}/({remove_letters(m.group(2))})"),
        
        # (a)/b
        (r'\(([A-Za-z0-9]+)\)/([A-Za-z0-9]+)',
         lambda m: f"({remove_letters(m.group(1))})/({remove_letters(m.group(2))})"),
        
        # Simple a/b
        (r'([A-Za-z0-9]+)/([A-Za-z0-9]+)',
         lambda m: f"{remove_letters(m.group(1))}/{remove_letters(m.group(2))}")
    ]

    result = diplotype
    for pattern, replacement in patterns:
        result = re.sub(pattern, replacement, result)
    
    return result

def diplotype_combinations_brackets(diplotype_pattern: str) -> list[list[str]]:
    # Regular expressions to match various patterns
    patterns = [
        r'\((\d+),(\d+)\)/\((\d+),(\d+)\)',  # (a,b)/(c,d)
        r'(\d+)/\((\d+),(\d+)\)',            # a/(b,c)
        r'\((\d+),(\d+)\)/(\d+)',            # (a,b)/c
        r'\((\d+),(\d+)\)/\((\d+)\)',        # (a,b)/(c)
        r'(\d+)/\((\d+)\)',                  # a/(b)
        r'\((\d+)\)/(\d+)',                  # (a)/b
        r'(\d+)/(\d+)'                       # a/b
    ]

    for pattern in patterns:
        match = re.match(pattern, diplotype_pattern)
        if match:
            numbers = list(map(int, match.groups()))
            
            if len(numbers) == 4:  # (a,b)/(c,d) case
                return list(map(list, product([str(numbers[0]), str(numbers[1])], [str(numbers[2]), str(numbers[3])])))
            elif len(numbers) == 3:
                if '/' in diplotype_pattern.split(')')[0]:  # Check if parentheses are on the right side
                    return [[str(numbers[0]), str(numbers[1])], [str(numbers[0]), str(numbers[2])]]
                else:
                    return [[str(numbers[0]), str(numbers[2])], [str(numbers[1]), str(numbers[2])]]
            elif len(numbers) == 2:
                return [[str(numbers[0]), str(numbers[1])]]

    raise ValueError(f"No matching pattern found for input: {diplotype_pattern}")

def categorize_alleles(allele: str) -> str:
    ### function to categorize alleles in sub and core alleles. 
    if "." in allele:
        return "sub-allele"
    if not "A" in allele.split("*")[1] and allele.split("*")[1][-1].isalpha():
        return "sub-allele"       
    return "core-allele"

def compare_versions_counts(df: pd.DataFrame) -> pd.DataFrame:
    versions = sorted([col for col in df.columns if col != 'Allele'])
    changes = []

    for i in range(len(versions) - 1):
        version_x = versions[i]
        version_y = versions[i + 1]
        
        alleles_x = set(df[df[version_x] == 1]['Allele'])
        alleles_y = set(df[df[version_y] == 1]['Allele'])
        
        removed = -len(alleles_x - alleles_y)
        added = len(alleles_y - alleles_x)
        
        source_value = df[df[version_x] == 1].shape[0]
        target_value = df[df[version_y] == 1].shape[0]
        
        if removed != 0:
            changes.append([version_x, version_y, 'Removed', removed, source_value, target_value])
        changes.append([version_x, version_y, 'Added', added, source_value, target_value])

    return pd.DataFrame(changes, columns=['source', 'target', 'change_type', 'value', 'source_value', 'target_value'])

def generate_sankey_data(data):
    nodes = list(pd.unique(data[['source', 'target']].values.ravel('K')))
    dummy_nodes = [f"Removed_{node}" for node in nodes]
    dummy_added = [f"Added_{node}" for node in nodes]
    all_nodes = nodes + dummy_nodes + dummy_added
    node_indices = {node: index for index, node in enumerate(all_nodes)}


    links = {
        "source": [],
        "target": [],
        "value": [],
        "label": [],
        "color": []
    }

    for i in range(len(nodes) - 1):
        source = nodes[i]
        target = nodes[i+1]

        dummy_target = f"Removed_{source}"
        dummy_added_source = f"Added_{target}"

        
        subset = data[(data['source'] == source) & (data['target'] == target)]

        
        if not subset.empty:
            source_value = subset['source_value'].iloc[0]
            target_value = subset['target_value'].iloc[0]
            added = subset[subset['change_type'] == 'Added']['value'].sum()
            removed = abs(subset[subset['change_type'] == 'Removed']['value']).sum()
            kept = source_value - removed

            # Kept alleles (first)
            # Added alleles (top)
            if added > 0:
                links["source"].append(node_indices[dummy_added_source])
                links["target"].append(node_indices[target])
                links["value"].append(added)
                links["label"].append("Added")
                links["color"].append("yellow")

            # Kept alleles (middle)
            if kept > 0:
                links["source"].append(node_indices[source])
                links["target"].append(node_indices[target])
                links["value"].append(kept)
                links["label"].append("Kept")
                links["color"].append("green")
            
            # Removed alleles (bottom)
            if removed > 0:
                links["source"].append(node_indices[source])
                links["target"].append(node_indices[dummy_target])
                links["value"].append(removed)
                links["label"].append("Removed")
                links["color"].append("red")
        else:
            # If no data, add a link with the previous total amount
            prev_total = data[data['source'] == source]['source_value'].iloc[0]
            links["source"].append(node_indices[source])
            links["target"].append(node_indices[target])
            links["value"].append(prev_total)
            links["label"].append("No Change")
            links["color"].append("gray")
        # Calculate positions
        n_real_nodes = len(nodes) // 2
        x_real = [i/(n_real_nodes-1) for i in range(n_real_nodes)]
        y_real = [0.3] * n_real_nodes

        x_dummy = x_real.copy()
        y_dummy = [0.7] * n_real_nodes

        # Combine positions
        x_positions = x_real + x_dummy
        y_positions = y_real + y_dummy

    return all_nodes, links, x_positions, y_positions

def create_sankey_plot(nodes, links, x_positions, y_positions):
    fig = go.Figure(data=[go.Sankey(
      node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = [label if not (label.startswith("Removed") or label.startswith("Added")) else "" for label in nodes],
      color = ["green" if not (n.startswith("Removed") or n.startswith("Added")) else "yellow" if n.startswith("Added") else "red" for n in nodes],
      x = x_positions,
      y = y_positions
    ),
    link = links)
    ])
    fig.update_layout(title_text="Major release Sankey Diagram", font_size=25)
    return fig


def add_star_prefix_to_list_elements(haplotype_lists: list[list[str]]):
    return [["*" + haplotype for haplotype in haplotype_list] for haplotype_list in haplotype_lists]

def create_list_of_haplotypes_from_diplotype(diplotype: str) -> list[list[str]]:
    if ("(" or ")") in diplotype:
        if ";" in diplotype:
            return process_diplotypes(diplotype, "both")
        else:
            return process_diplotypes(diplotype, "brackets")
    elif ";" in diplotype:
        return process_diplotypes(diplotype, "semicolon")
    
    return add_star_prefix_to_list_elements([diplotype.split("/")])
    
def process_diplotypes(diplotype: str, case: str) -> list[list[str]]:
    match case:
        case "brackets":
            return add_star_prefix_to_list_elements(diplotype_combinations_brackets(diplotype))

        case "semicolon":
            haplotype_list = diplotype.split(";")
            haplotype_list1 = haplotype_list[0].split("/")
            haplotype_list2 = haplotype_list[1].split("/")
            haplotype_lists = [haplotype_list1, haplotype_list2]
            return add_star_prefix_to_list_elements(haplotype_lists)

        case "both":
            haplotype_list = diplotype.split(";")
            result = []
            for i in range(len(haplotype_list)):
                if ("(" or ")") in haplotype_list[i]:
                    tmp_list = process_diplotypes(haplotype_list[i], "brackets")
                else:
                    tmp_list_2 = add_star_prefix_to_list_elements([haplotype_list[i].split("/")])

            if tmp_list is not None:
                result.extend(tmp_list)
            if tmp_list_2 is not None:
                result.extend(tmp_list_2)
            return result
        # Default case
        case _:
            raise ValueError(f"Unsupported case: {case}")

def predict_phenotypes(gene: str, diplotypes: list[list[str]]) -> list[str]:
    return [pypgx.predict_phenotype(gene, diplotype[0], diplotype[1]) for diplotype in diplotypes]
