import os
import requests
from Bio.Seq import Seq
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio.SeqUtils.MeltingTemp import Tm_NN

class ASOGenerator:
    # def __init__(self, transcript_id, aso_type="SSO", region="exon",
    #              target_exon=None, window_min=16, window_max=25,
    #              GC_min=40, GC_max=65, homopolymer_len=4,
    #              splice_window=5, gapmer_wing=5, gapmer_gap_min=7):
    #     self.transcript_id = transcript_id
    #     self.aso_type = aso_type
    #     self.region = region
    #     self.target_exon = target_exon
    #     self.window_min = window_min
    #     self.window_max = window_max
    #     self.GC_min = GC_min
    #     self.GC_max = GC_max
    #     self.homopolymer_len = homopolymer_len
    #     self.splice_window = splice_window
    #     self.gapmer_wing = gapmer_wing
    #     self.gapmer_gap_min = gapmer_gap_min

    #     self.sequence = ""
    #     self.exons = {}

    #     self.output_dir = os.path.expanduser("~/Downloads")
    #     os.makedirs(self.output_dir, exist_ok=True)

    # def __init__(self, transcript_id, aso_type="SSO", region="exon",
    #          target_exon=None, window_min=16, window_max=25,
    #          GC_min=40, GC_max=65, homopolymer_len=4,
    #          splice_window=5, gapmer_wing=5, gapmer_gap_min=7,
    #          genome_mode="local", output_dir="/tmp"):
    #     self.transcript_id = transcript_id
    #     self.aso_type = aso_type
    #     self.region = region
    #     self.target_exon = target_exon
    #     self.window_min = window_min
    #     self.window_max = window_max
    #     self.GC_min = GC_min
    #     self.GC_max = GC_max
    #     self.homopolymer_len = homopolymer_len
    #     self.splice_window = splice_window
    #     self.gapmer_wing = gapmer_wing
    #     self.gapmer_gap_min = gapmer_gap_min
    #     self.genome_mode = genome_mode

    #     self.sequence = ""
    #     self.exons = {}

    #     self.output_dir = output_dir
    #     os.makedirs(self.output_dir, exist_ok=True)

    def __init__(self, transcript_id=None, gene_name=None, ensembl_gene_id=None, 
                aso_type="SSO", region="exon", target_exon=None,
                window_min=16, window_max=25,
                GC_min=40, GC_max=65, homopolymer_len=4,
                splice_window=5, gapmer_wing=5, gapmer_gap_min=7,
                genome_mode="local", output_dir="/tmp"):
        
        self.transcript_id = transcript_id
        self.gene_name = gene_name
        self.ensembl_gene_id = ensembl_gene_id
        self.aso_type = aso_type
        self.region = region
        self.target_exon = target_exon
        self.window_min = window_min
        self.window_max = window_max
        self.GC_min = GC_min
        self.GC_max = GC_max
        self.homopolymer_len = homopolymer_len
        self.splice_window = splice_window
        self.gapmer_wing = gapmer_wing
        self.gapmer_gap_min = gapmer_gap_min
        self.genome_mode = genome_mode
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

        self.sequence = ""
        self.exons = {}
        self.chromosome = ""
        self.filter_counts = {}

        # Loading gene name from transcript
        if transcript_id and not gene_name:
            self.load_gene_name_from_transcript(transcript_id)
        elif gene_name and not transcript_id:
            self.load_mane_transcript(gene_name)

    def load_sequence(self):
        url = f"https://rest.ensembl.org/sequence/id/{self.transcript_id}?type=cdna"
        response = requests.get(url, headers={"Content-Type": "text/plain"})
        if response.ok:
            self.sequence = response.text.strip().upper().replace("T", "U")
        else:
            raise Exception(f"Ошибка: {response.status_code} — {response.text}")

    def load_exons(self):
        url = f"https://rest.ensembl.org/lookup/id/{self.transcript_id}?expand=1"
        response = requests.get(url, headers={"Content-Type": "application/json"})
        if response.ok:
            data = response.json()
            self.chromosome = data.get("seq_region_name", "")
            exons_list = data["Exon"]
            exons_list.sort(key=lambda x: x["start"])
            self.exons = {}
            for idx, exon in enumerate(exons_list, start=1):
                self.exons[idx] = {
                    "start": exon["start"],
                    "end": exon["end"],
                    "length": exon["end"] - exon["start"] + 1,
                }
        else:
            raise Exception(f"Error: {response.status_code} — {response.text}")

    def download_genomic_fragment(self):
        url = f"https://rest.ensembl.org/sequence/id/{self.transcript_id}?type=genomic"
        response = requests.get(url, headers={"Content-Type": "text/plain"})
        if response.ok:
            seq = response.text.strip()
            genome_path = os.path.join(self.output_dir, "genome.fa")
            with open(genome_path, "w") as f:
                f.write(f">genome_fragment\n{seq}")
            self.genome_fa = genome_path
        else:
            raise Exception(f"Genome loading error: {response.status_code}")

    # def build_bowtie2_index(self):
    #     index_prefix = os.path.join(self.output_dir, "genome_index")
    #     cmd = f"bowtie2-build {self.genome_fa} {index_prefix}"
    #     subprocess.run(cmd, shell=True, check=True)
    #     self.index_prefix = index_prefix

    def build_bowtie2_index(self):
        index_prefix = os.path.join(self.output_dir, "genome_index")
        cmd = f"bowtie2-build {self.genome_fa} {index_prefix}"
        subprocess.run(cmd, shell=True, check=True)
        self.index_prefix = index_prefix

    def target_region(self):
        if self.region == "full":
            return self.sequence
        elif self.region == "exon" and self.target_exon is not None:
            exon = self.exons[self.target_exon]
            exon_start_idx = sum(e["length"] for idx, e in self.exons.items() if idx < self.target_exon)
            exon_end_idx = exon_start_idx + exon["length"]
            return self.sequence[exon_start_idx:exon_end_idx]
        else:
            raise Exception("Invalid setting for region or target_exon")

    # def enumerate_windows(self, region_seq):
    #     fragments = []
    #     for w in range(self.window_min, self.window_max + 1):
    #         for i in range(len(region_seq) - w + 1):
    #             frag = region_seq[i:i + w]
    #             gc = self.gc_content(frag)
    #             if self.GC_min <= gc <= self.GC_max and not self.has_homopolymer(frag):
    #                 fragments.append({
    #                     "seq": frag,
    #                     "GC%": gc,
    #                     "length": len(frag),
    #                     "start": i,
    #                     "end": i + w
    #                 })
    #     return fragments

    def enumerate_windows(self, region_seq):
        total = 0
        passed_gc = 0
        passed_homopolymer = 0
        fragments = []

        for w in range(self.window_min, self.window_max + 1):
            for i in range(len(region_seq) - w + 1):
                total += 1
                frag = region_seq[i:i + w]
                gc = self.gc_content(frag)
                if self.GC_min <= gc <= self.GC_max:
                    passed_gc += 1
                    if not self.has_homopolymer(frag):
                        passed_homopolymer += 1
                        fragments.append({
                            "seq": frag,
                            "GC%": gc,
                            "length": len(frag),
                            "start": i,
                            "end": i + w
                        })
        
        self.filter_counts = {
            "Total": total,
            "GC filter": passed_gc,
            "Homopolymer filter": passed_homopolymer
        }
        return fragments

    def gc_content(self, seq):
        return (seq.count('G') + seq.count('C')) / len(seq) * 100

    def has_homopolymer(self, seq):
        return any(base * self.homopolymer_len in seq for base in "AUGC")

    def generate_aso(self, seq):
        temp_dna = seq.replace("U", "T")
        rc = str(Seq(temp_dna).reverse_complement())
        return rc.replace("T", "U")

    def write_fasta(self, asos):
        fasta_path = os.path.join(self.output_dir, "asos.fasta")
        with open(fasta_path, "w") as f:
            for idx, a in enumerate(asos):
                f.write(f">ASO_{idx}\n{self.generate_aso(a['seq']).replace('U','T')}\n")
        print(f"FASTA with ASO saved: {fasta_path}")
        self.fasta_path = fasta_path

    # def run_bowtie2(self):
    #     sam_path = os.path.join(self.output_dir, "bowtie2_out.sam")
    #     cmd = f"bowtie2 -f -x {self.index_prefix} -U {self.fasta_path} -S {sam_path} -p 4"
    #     subprocess.run(cmd, shell=True, check=True)
    #     self.sam_path = sam_path

    def run_bowtie2(self):
        if not hasattr(self, "index_prefix"):
            raise Exception("Bowtie2 index path is not set. Did you run build_bowtie2_index()?")

        sam_path = os.path.join(self.output_dir, "bowtie2_out.sam")
        cmd = f"bowtie2 -f -x {self.index_prefix} -U {self.fasta_path} -S {sam_path} -p 4 -k 10"
        subprocess.run(cmd, shell=True, check=True)
        self.sam_path = sam_path

    def parse_bowtie2_output(self):
        hits = {}
        with open(self.sam_path) as f:
            for line in f:
                if line.startswith('@'):
                    continue
                read_id = line.split('\t')[0]
                if read_id not in hits:
                    hits[read_id] = 0
                if line.split('\t')[2] != '*':
                    hits[read_id] += 1
        return hits

    # def run(self):
    #     self.load_sequence()
    #     self.load_exons()
    #     self.download_genomic_fragment()
    #     self.build_bowtie2_index()
    #     region_seq = self.target_region()
    #     print(f"Длина целевого региона: {len(region_seq)} nt")

    #     asos = self.enumerate_windows(region_seq)
    #     print(f"Найдено окон после фильтрации: {len(asos)}")

    #     self.write_fasta(asos)
    #     self.run_bowtie2()
    #     bowtie_hits = self.parse_bowtie2_output()

    #     # Формируем итоговую таблицу красиво
    #     records = []
    #     for idx, a in enumerate(asos):
    #         record = {
    #             "Original_seq": a["seq"],
    #             "ASO": self.generate_aso(a["seq"]),
    #             "GC%": a["GC%"],
    #             "Length": a["length"],
    #             "Start": a["start"],
    #             "End": a["end"],
    #             "Bowtie2_hits": bowtie_hits.get(f"ASO_{idx}", 0)
    #         }
    #         records.append(record)

    #     df = pd.DataFrame(records)
    #     csv_path = os.path.join(self.output_dir, "ASO_candidates.csv")
    #     df.to_csv(csv_path, index=False)
    #     print(f"Таблица готова: {csv_path}")
    # def calc_tm(self, seq):
    #     """Расчёт температуры плавления (Tm) по простой формуле Wallace Rule"""
    #     seq_dna = seq.replace("U", "T")
    #     return 2 * (seq_dna.count("A") + seq_dna.count("T")) + 4 * (seq_dna.count("G") + seq_dna.count("C"))
    
    def calc_tm(self, seq):
        """Calculate melting temperature using Nearest-Neighbor model (Biopython)."""
        seq_dna = seq.replace("U", "T")
        return Tm_NN(seq_dna)
    
    def run(self):
        self.load_sequence()
        self.load_exons()
        if self.target_exon and self.target_exon not in self.exons:
            raise ValueError(f"Exon {self.target_exon} not found in transcript {self.transcript_id}. Available exons: {list(self.exons.keys())}")

        if self.genome_mode == "local":
            self.download_genomic_fragment()
            self.build_bowtie2_index()
        elif self.genome_mode == "full":
            self.build_bowtie2_index()
        else:
            raise ValueError("genome_mode должен быть 'local' или 'full'")

        region_seq = self.target_region()
        print(f"Length of the target region: {len(region_seq)} nt")

        asos = self.enumerate_windows(region_seq)
        print(f"Number of windows found after filtering: {len(asos)}")

        self.write_fasta(asos)
        self.run_bowtie2()
        bowtie_hits = self.parse_bowtie2_output()

        # Формируем итоговую таблицу красиво
        records = []
        for idx, a in enumerate(asos):
            # record = {
            #     "Original_seq": a["seq"],
            #     "ASO": self.generate_aso(a["seq"]),
            #     "GC%": a["GC%"],
            #     "Length": a["length"],
            #     "Start": a["start"],
            #     "End": a["end"],
            #     "Tm": self.calc_tm(a["seq"]),
            #     "Bowtie2_hits": bowtie_hits.get(f"ASO_{idx}", 0)
            # }
            record = {
                "Gene": self.gene_name,
                "Transcript_ID": self.transcript_id,
                "Original_seq": a["seq"],
                "ASO": self.generate_aso(a["seq"]),
                "GC%": a["GC%"],
                "Length": a["length"],
                "Start": a["start"],
                "End": a["end"],
                "Tm": self.calc_tm(a["seq"]),
                "Bowtie2_hits": bowtie_hits.get(f"ASO_{idx}", 0),
                "ASO_type": self.aso_type,
                "Exon_number": self.target_exon,
                "Chromosome": self.chromosome,
                "Genome_start": self.exons[self.target_exon]["start"] + a["start"] if self.target_exon else None,
                "Genome_end": self.exons[self.target_exon]["start"] + a["end"] if self.target_exon else None
            }
            records.append(record)

        df = pd.DataFrame(records)
        passed_bowtie = df[df["Bowtie2_hits"] <= 1].shape[0]
        self.filter_counts["Bowtie2 filter"] = passed_bowtie
        csv_path = os.path.join(self.output_dir, "ASO_candidates.csv")
        df.to_csv(csv_path, index=False)
        print(f"Table is ready: {csv_path}")

    def plot_chromosome_map(self, df):
        # Plot exon length and ASO positions (relative coordinates within exon)
        exon = self.exons[self.target_exon]
        exon_len = exon["length"]

        fig, ax = plt.subplots(figsize=(10, 2))

        # Draw exon as grey block
        ax.add_patch(patches.Rectangle((0, 0.4), exon_len, 0.2, color='lightgrey'))
        ax.text(exon_len / 2, 0.7, f"Exon {self.target_exon}", ha="center", fontsize=10)

        # Draw ASOs as red blocks
        for idx, row in df.iterrows():
            ax.add_patch(patches.Rectangle((row["Start"], 0.1), row["Length"], 0.2, color='red'))

        ax.set_xlim(0, exon_len)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Relative Position within Exon")
        ax.set_yticks([])
        ax.set_title("ASO Positions on Target Exon")
        return fig
    
    def plot_filtering_barplot(self):
        """Plot of the number of candidates at each filtering stage"""
        fig, ax = plt.subplots()
        labels = list(self.filter_counts.keys())
        values = list(self.filter_counts.values())

        ax.bar(labels, values, color="skyblue")
        ax.set_title("ASO Filtering Pipeline")
        ax.set_ylabel("Number of Candidates")
        ax.set_xlabel("Filtering Step")
        plt.xticks(rotation=30)
        return fig
    
    def load_gene_name_from_transcript(self, transcript_id):
        url = f"https://rest.ensembl.org/lookup/id/{transcript_id}?"
        headers = {"Content-Type": "application/json"}
        response = requests.get(url, headers=headers)
        if response.ok:
            data = response.json()
            self.gene_name = data.get("display_name", "")
            self.ensembl_gene_id = data.get("Parent", "")
        else:
            raise Exception(f"Error loading gene name: {response.status_code} - {response.text}")

    # def load_mane_transcript(self, gene_name):
    #     url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
    #     headers = {"Content-Type": "application/json"}
    #     response = requests.get(url, headers=headers)

    #     if response.ok:
    #         data = response.json()
    #         print("Available transcripts:")
    #         for t in data.get("Transcript", []):
    #             print(f"{t['id']} — {t.get('version')} — {t.get('biotype')}")
    #         self.transcript_id = data.get("Transcript", [])[0]["id"]
    #         self.ensembl_gene_id = data.get("id", "")
    #     else:
    #         raise Exception(f"Ошибка загрузки транскрипта: {response.status_code} - {response.text}")
    def load_mane_transcript(self, gene_name):
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
        headers = {"Content-Type": "application/json"}
        response = requests.get(url, headers=headers)

        if response.ok:
            data = response.json()
            self.ensembl_gene_id = data.get("id", "")
            transcripts = data.get("Transcript", [])
            self.transcript_choices = [
                {
                    "display": f"{t.get('display_name', '')} ({t['id']})",
                    "id": t["id"]
                }
                for t in transcripts if t.get("biotype") == "protein_coding"
            ]
        else:
            raise Exception(f"Error loading transcript: {response.status_code} - {response.text}")

if __name__ == "__main__":
    gen = ASOGenerator(
        transcript_id=transcript_id,
        aso_type="SSO", #or GAPMER
        region="exon",
        target_exon=3,
        window_min=18,
        window_max=25,
        GC_min=40,
        GC_max=65
    )
    gen.run()