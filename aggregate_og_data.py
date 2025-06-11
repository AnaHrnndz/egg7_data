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
            "e6-seq2best_description.tsv",
            "e6-seq2bestname.tsv",
            "e6-seq2bigg.tsv",
            "e6.5-seq2card.tsv",
            "e6-seq2cazy.tsv",
        #    "e6.seq2go.tsv",
            "e6-seq2goslim.tsv",
            "e6.5-seq2kegg.tsv",
            "e6.5-seqs2pdb.tsv",
            "e6.5-seq2pfam_info.tsv",
            "e6-seq2smart.tsv",
        ]

        for filename in annotation_files:
            file_path = os.path.join(self.seq_annotation_dir, filename)
            key = filename.split("2")[1].replace('.tsv', '')
            print(key)
            
            if key == 'kegg':
                #kegg_annotations = ['knumber', 'reaction', 'pathway', 'module', 'enzyme', 'gsymbol', 'gname', 'desc']
                kegg_annotations = ['knumber', 'pathway', 'module', 'enzyme', 'brite', 'gsymbol', 'gname']

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

            # elif key == 'go':
                # self.annotations[key] = defaultdict(list)
                # self.annotations = parser_go(file_path, self.annotations)

            elif key == 'goslim':
                self.annotations[key] = defaultdict(list)
                self.annotations = parser_goslim(file_path, self.annotations)

            elif key == 'pdb':
                self.annotations[key] = defaultdict(list)
                self.annotations = parser_pdb(file_path, self.annotations)

            elif key == 'pfam_info':
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

        print(self.annotations.keys())
    
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
                    print(headers)
                    continue
                parts = line.strip().split("\t")
                og_entry = dict(zip(headers, parts))
                # Parse sequences into a list
                og_entry["Seqs"] = og_entry["Seqs"].split(",")
                self.ogs_data.append(og_entry)
        print("OGs data loaded.")

    def clean_dict(self, diccionario):
        """
        Elimina las claves cuyo valor esté vacío: None, '', [], {}, set()
        """
        
        clean_dict = defaultdict(dict)

        for k, val in diccionario.items():
            if k == 'fprof_sum' and isinstance(val, dict):
                clean_dict[k] = {f: v for f, v in val.items() if v}

            else:
                clean_dict[k] = val

        return clean_dict



        #return {k: v for k, v in diccionario['fprof_sum'].items() if v not in (None, '', [], {}, set())}


    def cross_link_information(self):
        """
        Cross-link OG information with sequence annotations.
        """
        print("Cross-linking OGs with sequence annotations...")
        ogs_with_annotations = []

        for og in self.ogs_data:
            
            seq_ids = og["Seqs"]
        
            spcies = set(s.split('.')[0] for s in seq_ids)
                
            clust_name = og["#OG_name"].split('@')[0]

            # Calculate top terms
            kegg_top_terms = self.calculate_top_terms(seq_ids, "knumber")
            kpathway_top_terms = self.calculate_top_terms(seq_ids, "pathway")
            kmodule_top_terms = self.calculate_top_terms(seq_ids, "module")
            kgsymbol_top_terms = self.calculate_top_terms(seq_ids, "gsymbol")
            kgname_top_terms = self.calculate_top_terms(seq_ids, "gname")
            best_name_top_term = self.calculate_top_terms(seq_ids, "bestname")
            bigg_top_term = self.calculate_top_terms(seq_ids, "bigg")
            cazy_top_term = self.calculate_top_terms(seq_ids, "cazy")
            card_top_term = self.calculate_top_terms(seq_ids, "card")
            #go_top_terms = self.calculate_top_terms(seq_ids, "go")
            goslim_top_term = self.calculate_top_terms(seq_ids, "goslim")
            pdb_top_term = self.calculate_top_terms(seq_ids, "pdb")
            pfam_top_term = self.calculate_top_terms(seq_ids, "pfam_arc")
            smart_top_term = self.calculate_top_terms(seq_ids, "smart_arc")

  
            if og["Inparalogs_Rate"] == '-':
                in_rate = '-'
            else: 
                in_rate = float(og["Inparalogs_Rate"])

            if og["SP_overlap_dup"] == '-':
                so_ovlap = '-'
            else:
                so_ovlap = float(og["SP_overlap_dup"])

            if og["OG_down"] == '-':
                nogdown = 0
            else:
                nogdown = len(og["OG_down"].split(','))

            if og["OG_up"] == '-':
                nogup = 0
            else:
                nogup = len(og["OG_up"].split(','))

            og_info = {
                "clust_name": clust_name,
                "n": og["#OG_name"],       # OG name
                "l": og["TaxoLevel"],      # Taxonomic level
                "sn": og["SciName_TaxoLevel"],  # Scientific name
                "an": og["AssocNode"],      # Associated node
                "ns": int(og["NumSP"]),     # Number of species
                "ch": og["OG_down"],        # OG down
                "par": og["OG_up"],          # OG up
                "npar": nogup,
                "nch": nogdown,
                "nm": int(og["NumSeqs"]), # Number of sequences
                "nrec": int(og["NumRecoverySeqs"]),  # Number of recovery sequences
                "ld": og["Lca_Dup"],        # LCA duplication
                "so": og["Species_Outliers"],  # Species outliers
                "nso": int(og["Num_SP_Outliers"]),  # Number of species outliers
                "ir": in_rate,  # Inparalogs rate
                "so_dup": so_ovlap,  # Species overlap duplication
                "seqs": seq_ids,       #list with annotations that will be added later         
                "spcs": list(spcies), 
                "rec_seqs": og["RecoverySeqs"],  # Recovery sequences
                "fprof_sum": {
                    "kegg_ko": kegg_top_terms,
                    "kegg_pathway": kpathway_top_terms,
                    "kegg_module": kmodule_top_terms,
                    "kegg_gsymbol": kgsymbol_top_terms,
                    "kegg_gname": kgname_top_terms,
                    #"go" : go_top_terms,
                    "best_name" : best_name_top_term,
                    "bigg" : bigg_top_term,
                    "cazy" : cazy_top_term,
                    "card" : card_top_term,
                    "GOslim" : goslim_top_term,
                    "pdb" : pdb_top_term,
                    "pfam_arch" : pfam_top_term,
                    "smart_arch" : smart_top_term
                }
            }

            # for sid in seq_ids:
                # seq_annotations = {key: self.annotations[key].get(sid, []) for key in self.annotations.keys()}
                # seq_info = {"id": sid, "ann": seq_annotations}
                # og_info["seqs"].append(seq_info)
            
            # for sid in seq_ids:
                # seq_annotations = {}
                # for key in self.annotations:
                    # annot = self.annotations[key].get(sid, [])
                    # if annot:
                        # seq_annotations[key] = annot
                
                # seq_info = {"id": sid, "ann": seq_annotations}
                # og_info["seqs"].append(seq_info)


            clean_og_info = self.clean_dict(og_info)
            ogs_with_annotations.append(clean_og_info)

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
            return list()

        term_counter = Counter()
        for seq_id in seq_ids:
            annotations = list(set(self.annotations[annotation_type].get(seq_id, [None])))
            term_counter.update(annotations)
        del term_counter[None]
        total = len(seq_ids)
        if total == 0:
            return list()

        top_terms = list()
        
        for term, count in term_counter.most_common(3):
            percentage = (count / total) * 100
            top_terms.append([term, f"{percentage:.2f}%"])
    
        
        return top_terms


if __name__ == "__main__":
    # Paths
    seq_annotation_dir = "/home/plaza/projects/eggnog6.5/analysis/seq_annotations"  # Directory containing annotation files
    ogs_file = os.path.join("/home/plaza/projects/eggnog6.5/analysis/og_delineation/", "e6.5-ogs_combined.tsv")  # Path to OGs file
    output_file = "/home/plaza/projects/eggnog6.5/analysis/ogs_func_annot/e6.5-fannot_ogs.jsonl"  # Output file

    # Preprocessor
    preprocessor = OGPreprocessor(seq_annotation_dir, ogs_file, output_file)
    preprocessor.run()
