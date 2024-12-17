import os
import json
import sys
from collections import defaultdict, Counter
from parser_data import parser_kegg, parser_best_desc, parser_best_name, parser_bigg, \
    parser_card, parser_cazy, parser_go, parser_goslim, parser_pdb, parser_pfam, parser_smart


class OGPreprocessor:
    def __init__(self, seq_annotation_dir: str, ogs_file: str, output_file: str):
        self.seq_annotation_dir = seq_annotation_dir
        self.ogs_file = ogs_file
        self.output_file = output_file
        self.annotations = {}
        self.ogs_data = []

    def load_annotations(self):
        """
        Load all sequence annotations from files in the annotation directory.
        """
        print("Loading sequence annotations...")
        annotation_files = [
            "e6.seq2best_description.tsv",
            "e6.seq2bestname.tsv",
            "e6.seq2bigg.tsv",
            "e6.seq2card.tsv",
            "e6.seq2cazy.tsv",
            "e6.seq2go.tsv",
            "e6.seq2goslim.tsv",
            "e6.seq2kegg.tsv",
            "e6.seq2pdb.tsv",
        #    "e6.seq2pfam.tsv",
            "e6.seq2pfam.no_clans_overlap.tsv",
            "e6.seq2smart.tsv",
        ]

        for filename in annotation_files:
            file_path = os.path.join(self.seq_annotation_dir, filename)
            key = filename.split(".")[1].split('2')[1] 
            print(key)
            
            if key == 'kegg':
                kegg_annotations = ['knumber', 'reaction', 'pathway', 'module', 'enzyme', 'gsymbol', 'gname', 'desc']

                for k in kegg_annotations:
                    self.annotations[k] = defaultdict(list)
                
                self.annotations = parser_kegg(file_path, self.annotations) 
                
            elif key == 'best_description':
                self.annotations[key] = defaultdict(list)
                self.annotations = parser_best_desc(file_path, self.annotations)
                
            elif key == 'bestname':
                self.annotations[key] = defaultdict(list)
                self.annotations = parser_best_name(file_path, self.annotations)
                
            elif key == 'bigg':
                self.annotations[key] = defaultdict(list)
                self.annotations = parser_bigg(file_path, self.annotations)

            elif key == 'card':
                self.annotations[key] = defaultdict(list)
                self.annotations = parser_card(file_path, self.annotations)

            elif key == 'cazy':
                self.annotations[key] = defaultdict(list)
                self.annotations = parser_cazy(file_path, self.annotations)

            elif key == 'go':
                self.annotations[key] = defaultdict(list)
                self.annotations = parser_go(file_path, self.annotations)

            elif key == 'goslim':
                self.annotations[key] = defaultdict(list)
                self.annotations = parser_goslim(file_path, self.annotations)

            elif key == 'pdb':
                self.annotations[key] = defaultdict(list)
                self.annotations = parser_pdb(file_path, self.annotations)

            elif key == 'pfam':
                self.annotations['pfam'] = defaultdict(list)
                self.annotations['pfam_arc'] = defaultdict(list)
                self.annotations = parser_pfam(file_path, self.annotations)

            elif key == 'smart':
                self.annotations['smart'] = defaultdict(list)
                self.annotations['smart_arc'] = defaultdict(list)
                self.annotations = parser_smart(file_path, self.annotations)
            

            else:
                self.annotations[key] = defaultdict(list)
                with open(file_path, "r") as file:
                    for line in file:
                        parts = line.strip().split("\t")
                        seq_id = parts[0]
                        values = parts[1:]
                        self.annotations[key][seq_id].extend(values)

        print("Annotations loaded.")



    def load_ogs(self):
        """
        Load OGs information from the OGs.tsv file.
        """
        print("Loading OGs data...")
        with open(self.ogs_file, "r") as file:
            headers = None
            for line in file:
                if line.startswith("#"):
                    headers = line.strip().split("\t")
                    continue
                parts = line.strip().split("\t")
                og_entry = dict(zip(headers, parts))
                # Parse sequences into a list
                og_entry["Seqs"] = og_entry["Seqs"].split(",")
                self.ogs_data.append(og_entry)
        print("OGs data loaded.")



    def cross_link_information(self):
        """
        Cross-link OG information with sequence annotations.
        """
        print("Cross-linking OGs with sequence annotations...")
        ogs_with_annotations = []

        for og in self.ogs_data:
            seq_ids = og["Seqs"]
            
            clust_name = og["#OG_name"].split('@')[0]

            # Calculate top terms
            kegg_top_terms = self.calculate_top_terms(seq_ids, "knumber")
            go_top_terms = self.calculate_top_terms(seq_ids, "go")
            best_name_top_term = self.calculate_top_terms(seq_ids, "bestname")
            bigg_top_term = self.calculate_top_terms(seq_ids, "bigg")
            cazy_top_term = self.calculate_top_terms(seq_ids, "cazy")
            card_top_term = self.calculate_top_terms(seq_ids, "card")
            go_top_term = self.calculate_top_terms(seq_ids, "go")
            goslim_top_term = self.calculate_top_terms(seq_ids, "goslim")
            pdb_top_term = self.calculate_top_terms(seq_ids, "pdb")
            pfam_top_term = self.calculate_top_terms(seq_ids, "pfam_arc")
            smart_top_term = self.calculate_top_terms(seq_ids, "smart_arc")
  
            
            og_info = {
                "clust_name": clust_name,
                "og": og["#OG_name"],       # OG name
                "tl": og["TaxoLevel"],      # Taxonomic level
                "sn": og["SciName_TaxoLevel"],  # Scientific name
                "an": og["AssocNode"],      # Associated node
                "sp": int(og["NumSP"]),     # Number of species
                "od": og["OG_down"],        # OG down
                "ou": og["OG_up"],          # OG up
                "nseq": int(og["NumSeqs"]), # Number of sequences
                "nrec": int(og["NumRecoverySeqs"]),  # Number of recovery sequences
                "ld": og["Lca_Dup"],        # LCA duplication
                "so": og["Species_Outliers"],  # Species outliers
                "nso": int(og["Num_SP_Outliers"]),  # Number of species outliers
                "ir": float(og["Inparalogs_Rate"]),  # Inparalogs rate
                "so_dup": float(og["SP_overlap_dup"]),  # Species overlap duplication
                "seqs": [],                 # Sequences with annotations
                "rec_seqs": og["RecoverySeqs"],  # Recovery sequences
                "top_kegg": kegg_top_terms,
                "top_go" : go_top_terms,
                "top_best_name" : best_name_top_term,
                "top_bigg" : bigg_top_term,
                "top_cazy" : cazy_top_term,
                "top_card" : card_top_term,
                "top_go" : go_top_term,
                "top_goslim" : goslim_top_term,
                "top_pdb" : pdb_top_term,
                "top_pfam" : pfam_top_term,
                "top_smart" : smart_top_term
            }

            for seq_id in og["Seqs"]:
                seq_annotations = {key: self.annotations[key].get(seq_id, []) for key in self.annotations}
                seq_info = {"id": seq_id, "ann": seq_annotations}
                og_info["seqs"].append(seq_info)

            ogs_with_annotations.append(og_info)

        print("Cross-linking complete.")
        return ogs_with_annotations



    def save_to_file(self, ogs_with_annotations):
        """
        Save the cross-linked OG data as JSON documents in a text file.
        """
        print("Saving cross-linked OG data to file...")
        with open(self.output_file, "w") as file:
            for og in ogs_with_annotations:
                file.write(json.dumps(og) + "\n")
        print(f"Data saved to {self.output_file}.")



    def run(self):
        """
        Run the full preprocessing pipeline.
        """
        self.load_annotations()
        self.load_ogs()
        ogs_with_annotations = self.cross_link_information()
        self.save_to_file(ogs_with_annotations)



    def calculate_top_terms(self, seq_ids, annotation_type):
        if annotation_type not in self.annotations:
            print(f"Warning: Annotation type '{annotation_type}' not found.")
            return []

        term_counter = Counter()
        for seq_id in seq_ids:
            annotations = list(set(self.annotations[annotation_type].get(seq_id, [None])))
            term_counter.update(annotations)
        del term_counter[None]
        total = len(seq_ids)
        if total == 0:
            return []

        top_terms = []
        
        for term, count in term_counter.most_common(3):
            percentage = (count / total) * 100
            top_terms.append([term, f"{percentage:.2f}%"])
    
        
        return top_terms


if __name__ == "__main__":
    # Paths
    seq_annotation_dir = "/data/projects/Egg7_data/raw_annotations"  # Directory containing annotation files
    ogs_file = os.path.join("/data/projects/Egg7_data/raw_results", "OGs.tsv")  # Path to OGs file
    output_file = "/data/projects/Egg7_data/process_results/output_ogs_with_annotations.jsonl"  # Output file

    # Preprocessor
    preprocessor = OGPreprocessor(seq_annotation_dir, ogs_file, output_file)
    preprocessor.run()
