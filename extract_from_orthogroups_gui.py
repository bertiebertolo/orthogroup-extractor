#!/usr/bin/env python3
"""
GUI: extract_from_orthogroups_gui.py

Paths-only Tkinter GUI that:
 - reads an OrthoFinder Orthogroups.tsv
 - reads a CSV of queries (select column)
 - selects source species (where CSV IDs originate) and target species (or "All")
 - auto-detects whether CSV column contains protein IDs (with suffix) or gene-prefixes
 - MATCHES using exact ID or gene-prefix (auto-enabled when appropriate)
 - DOES NOT require an external id-map file (id-map support removed)
 - Writes outputs:
     * For specific target species:
         <out_prefix>_<target>_target.tsv  (query_id, orthogroup, <target>_members)
         <out_prefix>_<target>_not_found.tsv (queries not found)
     * If target == "All":
         <out_prefix>_all_wide.tsv  (single wide file with one column per species EXCEPT source species)
         <out_prefix>_all_not_found.tsv (queries not found)
 - Does not include any original CSV columns in outputs.
 - Does not produce the long TSV anymore.
"""
import tkinter as tk
from tkinter import ttk, messagebox
import tkinter.font as tkfont
import csv
import re
from collections import OrderedDict, defaultdict
import threading
import queue
import os

# optional copy support
try:
    import pyperclip
    def pyperclip_copy(x):
        pyperclip.copy(x)
except Exception:
    def pyperclip_copy(x):
        return

# ---------------------
# Helpers
# ---------------------
def normalize_protid_header(hdr):
    return hdr.split()[0] if hdr else ""

def gene_prefix_from_protid(protid):
    token = normalize_protid_header(protid)
    # include hyphen in split so "AAEL027978-PA" -> "AAEL027978"
    return re.split(r'[-|._:]', token)[0]

def sanitize_fname(s):
    return re.sub(r'\s+', '_', s)

def parse_orthogroups_full(path, log_fn=None):
    og_rows = OrderedDict()
    prot_to_og = {}
    prot_to_species = {}
    geneprefix_index = defaultdict(set)
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        if len(header) < 2:
            raise RuntimeError("Orthogroups.tsv header looks wrong.")
        species = header[1:]
        if log_fn:
            log_fn(f"Orthogroups header parsed: {len(species)} species.")
        for i, row in enumerate(reader, start=1):
            if not row:
                continue
            og = row[0]
            cells = row[1:]
            og_rows[og] = cells
            for sp_name, cell in zip(species, cells):
                if not cell:
                    continue
                ids = [s.strip() for s in re.split(r'[;,]', cell) if s.strip()]
                for pid in ids:
                    normalized = normalize_protid_header(pid)
                    prot_to_og[normalized] = og
                    prot_to_species[normalized] = sp_name
                    gp = gene_prefix_from_protid(normalized)
                    geneprefix_index[gp].add(normalized)
            if log_fn and (i % 5000 == 0):
                log_fn(f"Parsed {i} orthogroups...")
    return species, og_rows, prot_to_og, prot_to_species, geneprefix_index

def sniff_csv_delimiter(sample_bytes):
    try:
        sample = sample_bytes.decode('utf-8', errors='replace')
    except AttributeError:
        sample = str(sample_bytes)
    sniffer = csv.Sniffer()
    try:
        dialect = sniffer.sniff(sample)
        return dialect.delimiter
    except Exception:
        for d in [',', '\t', ';']:
            if d in sample:
                return d
        return ','

def read_csv_header_and_preview(path, n_preview=20):
    with open(path, 'rb') as fh:
        sample = fh.read(8192)
    delim = sniff_csv_delimiter(sample)
    with open(path, newline='', encoding='utf-8', errors='replace') as fh:
        reader = csv.reader(fh, delimiter=delim)
        try:
            header = next(reader)
        except StopIteration:
            header = []
        preview_rows = []
        for i, row in enumerate(reader):
            if i >= n_preview:
                break
            preview_rows.append(row)
    return header, preview_rows, delim

def read_csv_rows(path):
    header, _, delim = read_csv_header_and_preview(path, n_preview=1)
    rows = []
    with open(path, newline='', encoding='utf-8', errors='replace') as fh:
        reader = csv.reader(fh, delimiter=delim)
        try:
            header = next(reader)
        except StopIteration:
            return [], [], delim
        for row in reader:
            rows.append([c for c in row])
    return header, rows, delim

def extract_column_from_rows(rows, header, column_name_or_index):
    if isinstance(column_name_or_index, int):
        col_idx = column_name_or_index
    else:
        try:
            col_idx = header.index(column_name_or_index)
        except ValueError:
            try:
                col_idx = int(column_name_or_index)
            except Exception:
                raise ValueError("Column name not found in CSV header.")
    out = []
    for row in rows:
        val = ""
        if col_idx < len(row):
            val = row[col_idx].strip()
        out.append((row, val))
    return out

def detect_column_id_style_from_rows(rows_with_vals, sample_n=200):
    vals = [val for (_, val) in rows_with_vals[:sample_n] if val]
    if not vals:
        return {'protein_frac': 0.0, 'prefix_frac': 0.0, 'mostly_prefix': False, 'sample_count':0}
    prot_count = 0
    prefix_count = 0
    prot_re = re.compile(r'.+-[A-Za-z0-9]+$')
    prefix_re = re.compile(r'^[A-Za-z]{1,}\d{3,}[A-Za-z0-9]*$')
    for v in vals:
        if prot_re.search(v):
            prot_count += 1
        elif prefix_re.match(v):
            prefix_count += 1
    total = len(vals)
    prot_frac = prot_count / total
    prefix_frac = prefix_count / total
    mostly_prefix = (prefix_frac >= 0.6 and prot_frac < 0.4)
    return {'protein_frac': prot_frac, 'prefix_frac': prefix_frac, 'mostly_prefix': mostly_prefix, 'sample_count': total}

def parse_members_cell(cell):
    if not cell:
        return []
    return [s.strip() for s in re.split(r'[;,]', cell) if s.strip()]

def format_ids(ids):
    return "|".join(ids)

# ---------------------
# GUI
# ---------------------
class App:
    def __init__(self, root):
        self.root = root
        root.title("OrthoFinder CSV -> Target (no id-map, paths-only)")

        self.log_q = queue.Queue()

        frm = ttk.Frame(root, padding=10)
        frm.grid(sticky="nsew")
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)

        # layout: left inputs, right ID preview
        left = ttk.Frame(frm)
        left.grid(row=0, column=0, sticky="nsew")
        right = ttk.Frame(frm)
        right.grid(row=0, column=1, sticky="nsew", padx=(8,0))
        frm.columnconfigure(0, weight=3)
        frm.columnconfigure(1, weight=2)

        # Inputs (left)
        ttk.Label(left, text="Orthogroups.tsv (full path):").grid(row=0, column=0, sticky="w")
        self.ortho_entry = ttk.Entry(left, width=90)
        self.ortho_entry.grid(row=0, column=1, sticky="we", padx=5)
        ttk.Button(left, text="Load species", command=self.load_species).grid(row=0, column=2, padx=(6,0))

        ttk.Label(left, text="Source species (CSV IDs come from):").grid(row=1, column=0, sticky="w")
        self.source_species_var = tk.StringVar()
        self.source_species_combo = ttk.Combobox(left, textvariable=self.source_species_var, state="readonly", width=60)
        self.source_species_combo.grid(row=1, column=1, sticky="we", padx=5)

        ttk.Label(left, text="Target species (find orthologs in) or All:").grid(row=2, column=0, sticky="w")
        self.target_species_var = tk.StringVar()
        self.target_species_combo = ttk.Combobox(left, textvariable=self.target_species_var, state="readonly", width=60)
        self.target_species_combo.grid(row=2, column=1, sticky="we", padx=5)

        ttk.Label(left, text="CSV queries file (full path):").grid(row=3, column=0, sticky="w", pady=(8,0))
        self.csv_entry = ttk.Entry(left, width=90)
        self.csv_entry.grid(row=3, column=1, sticky="we", padx=5, pady=(8,0))
        ttk.Button(left, text="Load CSV columns", command=self.load_csv_columns).grid(row=3, column=2, padx=(6,0))

        ttk.Label(left, text="CSV column:").grid(row=4, column=0, sticky="w", pady=(6,0))
        self.csv_col_var = tk.StringVar()
        self.csv_col_combo = ttk.Combobox(left, textvariable=self.csv_col_var, state="readonly", width=60)
        self.csv_col_combo.grid(row=4, column=1, sticky="we", padx=5)
        self.csv_col_combo.bind("<<ComboboxSelected>>", self.on_column_selected)

        ttk.Label(left, text="CSV preview (first values):").grid(row=5, column=0, sticky="nw", pady=(6,0))
        self.csv_preview = tk.Text(left, width=80, height=6)
        self.csv_preview.grid(row=5, column=1, columnspan=2, sticky="we", padx=5, pady=(6,0))

        # columns-to-keep selector (multi-select)
        ttk.Label(left, text="Columns to keep in output (select, Ctrl+click):").grid(row=6, column=0, sticky="w", pady=(6,0))
        self.keep_cols_listbox = tk.Listbox(left, selectmode="multiple", height=6, exportselection=False)
        self.keep_cols_listbox.grid(row=6, column=1, columnspan=2, sticky="we", padx=5)

        # match-prefix option (auto-enabled) - using tk.Checkbutton with wrapping to avoid clipping
        bfont = tkfont.Font(family="Helvetica", size=11)
        self.match_prefix_var = tk.BooleanVar(value=False)
        self.match_prefix_cb = tk.Checkbutton(
            left,
            text="Match gene-prefix (auto-enabled if CSV looks like gene IDs)",
            variable=self.match_prefix_var,
            font=bfont,
            anchor="w",
            justify="left",
            wraplength=600
        )
        # allow it to span both columns so text can use full width
        self.match_prefix_cb.grid(row=7, column=1, columnspan=2, sticky="we", pady=(6,0))

        ttk.Label(left, text="Output folder (full path):").grid(row=8, column=0, sticky="w", pady=(6,0))
        self.out_entry = ttk.Entry(left, width=90)
        self.out_entry.grid(row=8, column=1, sticky="we", padx=5, pady=(6,0))

        ttk.Label(left, text="Output prefix:").grid(row=9, column=0, sticky="w", pady=(6,0))
        self.out_prefix_entry = ttk.Entry(left, width=40)
        self.out_prefix_entry.insert(0, "orthologs")
        self.out_prefix_entry.grid(row=9, column=1, sticky="w", padx=5, pady=(6,0))

        self.run_btn = ttk.Button(left, text="Run", command=self.start_run)
        self.run_btn.grid(row=10, column=1, pady=(12,0), sticky="w")

        ttk.Label(left, text="Log / progress:").grid(row=11, column=0, sticky="nw", pady=(8,0))
        self.log_text = tk.Text(left, width=100, height=12, wrap="none")
        self.log_text.grid(row=11, column=1, columnspan=2, sticky="nsew", padx=5, pady=(8,0))
        left.rowconfigure(11, weight=1)
        left.columnconfigure(1, weight=1)

        # Right column: ID preview / list
        ttk.Label(right, text="IDs to be analyzed (unique preview):").grid(row=0, column=0, sticky="w")
        self.id_count_var = tk.StringVar(value="0 unique IDs")
        ttk.Label(right, textvariable=self.id_count_var).grid(row=0, column=1, sticky="w", padx=(6,0))

        ttk.Label(right, text="Filter:").grid(row=1, column=0, sticky="w")
        self.filter_entry = ttk.Entry(right, width=30)
        self.filter_entry.grid(row=1, column=1, sticky="we", padx=(6,0))
        ttk.Button(right, text="Apply", command=self.apply_filter).grid(row=1, column=2, padx=(6,0))

        self.ids_listbox = tk.Listbox(right, height=30)
        self.ids_listbox.grid(row=2, column=0, columnspan=3, sticky="nsew", pady=(6,0))
        right.rowconfigure(2, weight=1)
        right.columnconfigure(1, weight=1)
        self.ids_listbox.bind("<Double-Button-1>", self.on_id_double_click)

        ttk.Button(right, text="Copy selected ID", command=self.copy_selected_id).grid(row=3, column=2, sticky="e", pady=(6,0))
        ttk.Button(right, text="Refresh list", command=self.refresh_id_list).grid(row=3, column=0, sticky="w", pady=(6,0))

        # bind resize so checkbox wraplength follows window size
        self.root.bind("<Configure>", self._update_checkbox_wrap)
        self.root.after(100, self._update_checkbox_wrap)

        self.root.after(200, self.poll_log)

        # internal storage
        self.species = []
        self.og_rows = None
        self.prot_to_og = None
        self.prot_to_species = None
        self.geneprefix_index = None

        self.csv_header = []
        self.csv_rows = []
        self.rows_with_vals = []   # list of (row, value)
        self.unique_ids = []       # in-order unique list
        self.filtered_ids = []

    # UI actions
    def load_species(self):
        path = self.ortho_entry.get().strip()
        if not path or not os.path.isfile(path):
            messagebox.showerror("Error", "Please enter a valid Orthogroups.tsv full file path first.")
            return
        try:
            self.species, self.og_rows, self.prot_to_og, self.prot_to_species, self.geneprefix_index = parse_orthogroups_full(path, log_fn=self.log)
            # include 'All' as an option in target
            vals = ["All"] + list(self.species)
            self.source_species_combo['values'] = list(self.species)
            self.target_species_combo['values'] = vals
            if self.species:
                self.source_species_combo.current(0)
                self.target_species_combo.current(0)  # default to All
            self.log(f"Loaded Orthogroups.tsv: {len(self.species)} species available.")
        except Exception as e:
            messagebox.showerror("Error parsing Orthogroups.tsv", str(e))

    def load_csv_columns(self):
        path = self.csv_entry.get().strip()
        if not path or not os.path.isfile(path):
            messagebox.showerror("Error", "Please enter a valid CSV full file path first.")
            return
        try:
            header, rows, delim = read_csv_rows(path)
            if not header:
                messagebox.showerror("Error", "CSV appears to be empty or malformed.")
                return
            self.csv_header = header
            self.csv_rows = rows
            self.csv_col_combo['values'] = header
            self.csv_col_combo.current(0)
            # populate keep-columns listbox (default: select all)
            self.keep_cols_listbox.delete(0, tk.END)
            for i, col in enumerate(header):
                self.keep_cols_listbox.insert(tk.END, col)
                self.keep_cols_listbox.selection_set(i)
            preview_rows = rows[:30]
            preview_text = f"Detected delimiter: '{delim}'\nHeader:\n{header}\n\nFirst {len(preview_rows)} rows:\n"
            for row in preview_rows:
                preview_text += str(row) + "\n"
            self.csv_preview.delete(1.0, tk.END)
            self.csv_preview.insert(tk.END, preview_text)
            self.log(f"Loaded CSV header with {len(header)} columns. Preview shown.")

            # Auto-detect ID style for selected column
            selected_col = self.csv_col_combo.get()
            if selected_col:
                rows_with_vals = extract_column_from_rows(rows, header, selected_col)
                self.rows_with_vals = rows_with_vals
                detect = detect_column_id_style_from_rows(rows_with_vals, sample_n=200)
                if detect.get('mostly_prefix'):
                    self.match_prefix_var.set(True)
                    self.log(f"Auto-detected CSV column '{selected_col}' as gene-prefix style (prefix_frac={detect['prefix_frac']:.2f}); enabled Match gene-prefix.")
                else:
                    self.log(f"CSV column detection: protein_frac={detect['protein_frac']:.2f}, prefix_frac={detect['prefix_frac']:.2f} (no auto-change).")

                # build unique id list and show it
                self.build_unique_ids()
                self.refresh_id_list()
        except Exception as e:
            messagebox.showerror("Error reading CSV", str(e))

    def on_column_selected(self, event=None):
        col = self.csv_col_combo.get()
        if not self.csv_rows or not self.csv_header:
            return
        try:
            self.rows_with_vals = extract_column_from_rows(self.csv_rows, self.csv_header, col)
            detect = detect_column_id_style_from_rows(self.rows_with_vals, sample_n=200)
            if detect.get('mostly_prefix'):
                self.match_prefix_var.set(True)
                self.log(f"Auto-detected CSV column '{col}' as gene-prefix style; enabled Match gene-prefix.")
            else:
                self.log(f"CSV column detection: protein_frac={detect['protein_frac']:.2f}, prefix_frac={detect['prefix_frac']:.2f}")

            self.build_unique_ids()
            self.refresh_id_list()
        except Exception as e:
            messagebox.showerror("Error selecting column", str(e))

    def build_unique_ids(self):
        seen = set()
        uniq = []
        for _, val in self.rows_with_vals:
            if not val:
                continue
            v = val.strip()
            if v not in seen:
                seen.add(v)
                uniq.append(v)
        self.unique_ids = uniq
        self.filtered_ids = list(self.unique_ids)
        self.id_count_var.set(f"{len(self.unique_ids)} unique IDs")

    def refresh_id_list(self):
        self.ids_listbox.delete(0, tk.END)
        for v in self.filtered_ids:
            self.ids_listbox.insert(tk.END, v)

    def apply_filter(self):
        txt = self.filter_entry.get().strip().lower()
        if not txt:
            self.filtered_ids = list(self.unique_ids)
        else:
            self.filtered_ids = [v for v in self.unique_ids if txt in v.lower()]
        self.refresh_id_list()
        self.id_count_var.set(f"{len(self.unique_ids)} unique IDs ({len(self.filtered_ids)} shown)")

    def on_id_double_click(self, event=None):
        sel = self.ids_listbox.curselection()
        if not sel:
            return
        v = self.ids_listbox.get(sel[0])
        try:
            pyperclip_copy(v)
            self.log(f"Copied to clipboard: {v}")
        except Exception:
            self.log("Copy to clipboard not available on this system.")

    def copy_selected_id(self):
        sel = self.ids_listbox.curselection()
        if not sel:
            messagebox.showinfo("Copy", "No ID selected.")
            return
        v = self.ids_listbox.get(sel[0])
        try:
            pyperclip_copy(v)
            messagebox.showinfo("Copy", f"Copied to clipboard: {v}")
        except Exception:
            messagebox.showinfo("Copy", "Copy not available on this system.")

    def get_selected_keep_columns(self):
        sel = self.keep_cols_listbox.curselection()
        if not sel:
            return []
        return [self.keep_cols_listbox.get(i) for i in sel]

    def _update_checkbox_wrap(self, event=None):
        # adjust wraplength so the checkbox label reflows with window size
        try:
            # set wrap to ~65% of the root width but not less than 200px
            wrap_px = max(200, int(self.root.winfo_width() * 0.65))
            self.match_prefix_cb.config(wraplength=wrap_px)
        except Exception:
            pass

    def log(self, msg):
        self.log_q.put(msg)

    def poll_log(self):
        try:
            while True:
                msg = self.log_q.get_nowait()
                self.log_text.insert(tk.END, msg + "\n")
                self.log_text.see(tk.END)
        except queue.Empty:
            pass
        self.root.after(200, self.poll_log)

    def start_run(self):
        ortho = self.ortho_entry.get().strip()
        csv_file = self.csv_entry.get().strip()
        csv_col = self.csv_col_combo.get().strip()
        source_species = self.source_species_var.get().strip()
        target_species = self.target_species_combo.get().strip()
        out_folder = self.out_entry.get().strip() or os.getcwd()
        out_prefix = self.out_prefix_entry.get().strip() or "orthologs"
        match_prefix = bool(self.match_prefix_var.get())
        keep_cols = self.get_selected_keep_columns()

        if not ortho or not os.path.isfile(ortho):
            messagebox.showerror("Input error", "Please enter a valid Orthogroups.tsv file path.")
            return
        if not csv_file or not os.path.isfile(csv_file):
            messagebox.showerror("Input error", "Please enter a valid CSV queries file path.")
            return
        if not csv_col:
            messagebox.showerror("Input error", "Please load CSV columns and select a column.")
            return
        if not source_species:
            messagebox.showerror("Input error", "Please select source species (Load species first).")
            return
        if not out_folder or not os.path.isdir(out_folder):
            messagebox.showerror("Input error", "Please enter a valid output folder path.")
            return

        # ensure orthogroups parsed
        if self.og_rows is None:
            try:
                self.species, self.og_rows, self.prot_to_og, self.prot_to_species, self.geneprefix_index = parse_orthogroups_full(ortho, log_fn=self.log)
            except Exception as e:
                messagebox.showerror("Error parsing Orthogroups.tsv", str(e))
                return

        self.run_btn.config(state="disabled")
        self.log("Starting ortholog extraction (no id-map)...")
        t = threading.Thread(target=self.worker_run, args=(csv_file, csv_col, source_species, target_species, match_prefix, out_folder, out_prefix, keep_cols))
        t.daemon = True
        t.start()

    def worker_run(self, csv_file, csv_col, source_species, target_species, match_prefix, out_folder, out_prefix, keep_cols):
        try:
            header, rows, delim = read_csv_rows(csv_file)
            if not header:
                self.log("CSV appears empty; aborting.")
                messagebox.showerror("Error", "CSV appears empty.")
                self.run_btn_after_enable()
                return

            rows_with_vals = extract_column_from_rows(rows, header, csv_col)
            queries = rows_with_vals  # list of (row, value)
            total_queries = len(queries)
            self.log(f"Loaded {total_queries} queries from CSV column '{csv_col}' (including blank rows).")

            target_is_all = (str(target_species).strip().lower() == "all")
            safe_out_prefix = sanitize_fname(out_prefix)

            # Determine which indices from header correspond to keep_cols (for appending rows)
            # OVERRIDE: do not include any input CSV columns in outputs
            keep_cols = []
            keep_indices = [header.index(c) for c in keep_cols if c in header]

            if target_is_all:
                out_fname = os.path.join(out_folder, f"{safe_out_prefix}_all_wide.tsv")
                not_found_fname = os.path.join(out_folder, f"{safe_out_prefix}_all_not_found.tsv")
                # build species list excluding source_species
                species_excluding_source = [s for s in self.species if s != source_species]
                with open(out_fname, "w", newline='') as wfh, open(not_found_fname, "w", newline='') as nfh:
                    writer = csv.writer(wfh, delimiter="\t")
                    n_writer = csv.writer(nfh, delimiter="\t")
                    wide_fieldnames = ["query_id", "orthogroup"] + species_excluding_source
                    writer.writerow(wide_fieldnames)
                    n_writer.writerow(["query_id", "reason"])

                    count_not_found = 0
                    count_matched = 0
                    count_matched_og = 0

                    for i, (row, q) in enumerate(queries, start=1):
                        if i % 200 == 0:
                            self.log(f"Processing {i}/{total_queries} queries...")
                        q = (q or "").strip()
                        matched_proteins = set()
                        # exact matches (raw q or normalized) restricted to source_species
                        if q and q in self.prot_to_og and self.prot_to_species.get(q) == source_species:
                            matched_proteins.add(q)
                        normq = normalize_protid_header(q) if q else ""
                        if normq and normq in self.prot_to_og and self.prot_to_species.get(normq) == source_species:
                            matched_proteins.add(normq)
                        if match_prefix and q:
                            candidates = self.geneprefix_index.get(q, set())
                            for pid in candidates:
                                if self.prot_to_species.get(pid) == source_species:
                                    matched_proteins.add(pid)
                            qtok = re.split(r'[-|._:]', normq)[0] if normq else ""
                            if qtok and qtok != q:
                                candidates2 = self.geneprefix_index.get(qtok, set())
                                for pid in candidates2:
                                    if self.prot_to_species.get(pid) == source_species:
                                        matched_proteins.add(pid)

                        if not matched_proteins:
                            count_not_found += 1
                            n_writer.writerow([q, "NOT_FOUND"])
                            writer.writerow([q, ""] + [""]*len(species_excluding_source))
                            continue

                        ogs = {self.prot_to_og[p] for p in matched_proteins if p in self.prot_to_og}
                        if not ogs:
                            count_not_found += 1
                            n_writer.writerow([q, "NOT_FOUND"])
                            writer.writerow([q, ""] + [""]*len(species_excluding_source))
                            continue

                        count_matched += 1
                        count_matched_og += len(ogs)

                        for og in sorted(ogs):
                            members_row = self.og_rows.get(og)
                            wide_cells = []
                            if members_row is not None:
                                for sp in species_excluding_source:
                                    try:
                                        idx_sp = self.species.index(sp)
                                        cell = members_row[idx_sp] if idx_sp < len(members_row) else ""
                                    except ValueError:
                                        cell = ""
                                    ids = parse_members_cell(cell)
                                    wide_cells.append(format_ids(ids))
                            else:
                                wide_cells = [""] * len(species_excluding_source)
                            writer.writerow([q, og] + wide_cells)

                    summary_lines = [
                        f"Total queries processed: {total_queries}",
                        f"Queries with any OG match in source species: {count_matched}",
                        f"Total OG-rows written (sum over matching OGs): {count_matched_og}",
                        f"Queries NOT_FOUND: {count_not_found}"
                    ]
                    for ln in summary_lines:
                        self.log(ln)
                    messagebox.showinfo("Run summary", "\n".join(summary_lines))
                    self.log(f"Done. Wrote:\n  {out_fname}\n  {not_found_fname}")
            else:
                # specific target species: output one file and a not-found file
                safe_target = sanitize_fname(target_species)
                prefix = f"{safe_out_prefix}_{safe_target}"
                target_out = os.path.join(out_folder, f"{prefix}_target.tsv")
                not_found_out = os.path.join(out_folder, f"{prefix}_not_found.tsv")

                with open(target_out, "w", newline='') as tfh, open(not_found_out, "w", newline='') as nfh:
                    t_writer = csv.writer(tfh, delimiter="\t")
                    n_writer = csv.writer(nfh, delimiter="\t")
                    t_writer.writerow(["query_id", "orthogroup", f"{target_species}_members"])
                    n_writer.writerow(["query_id", "reason"])

                    count_not_found = 0
                    count_matched = 0
                    count_matched_target_nonempty = 0
                    count_matched_og = 0

                    for i, (row, q) in enumerate(queries, start=1):
                        if i % 200 == 0:
                            self.log(f"Processing {i}/{total_queries} queries...")
                        q = (q or "").strip()
                        matched_proteins = set()
                        if q and q in self.prot_to_og and self.prot_to_species.get(q) == source_species:
                            matched_proteins.add(q)
                        normq = normalize_protid_header(q) if q else ""
                        if normq and normq in self.prot_to_og and self.prot_to_species.get(normq) == source_species:
                            matched_proteins.add(normq)
                        if match_prefix and q:
                            candidates = self.geneprefix_index.get(q, set())
                            for pid in candidates:
                                if self.prot_to_species.get(pid) == source_species:
                                    matched_proteins.add(pid)
                            qtok = re.split(r'[-|._:]', normq)[0] if normq else ""
                            if qtok and qtok != q:
                                candidates2 = self.geneprefix_index.get(qtok, set())
                                for pid in candidates2:
                                    if self.prot_to_species.get(pid) == source_species:
                                        matched_proteins.add(pid)

                        if not matched_proteins:
                            count_not_found += 1
                            n_writer.writerow([q, "NOT_FOUND"])
                            t_writer.writerow([q, "", "NOT_FOUND"])
                            continue

                        ogs = {self.prot_to_og[p] for p in matched_proteins if p in self.prot_to_og}
                        if not ogs:
                            count_not_found += 1
                            n_writer.writerow([q, "NOT_FOUND"])
                            t_writer.writerow([q, "", "NOT_FOUND"])
                            continue

                        count_matched += 1
                        count_matched_og += len(ogs)

                        for og in sorted(ogs):
                            members_row = self.og_rows.get(og)
                            if members_row is None:
                                target_cell = ""
                            else:
                                try:
                                    idx_sp = self.species.index(target_species)
                                    target_cell = members_row[idx_sp] if idx_sp < len(members_row) else ""
                                except ValueError:
                                    target_cell = ""
                            target_ids = parse_members_cell(target_cell)
                            target_str = format_ids(target_ids)

                            if target_str:
                                count_matched_target_nonempty += 1

                            t_writer.writerow([q, og, target_str])

                summary_lines = [
                    f"Total queries processed: {total_queries}",
                    f"Queries with any OG match in source species: {count_matched}",
                    f"Queries with non-empty members in target species: {count_matched_target_nonempty}",
                    f"Total OG-rows written (sum over matching OGs): {count_matched_og}",
                    f"Queries NOT_FOUND: {count_not_found}"
                ]
                for ln in summary_lines:
                    self.log(ln)
                messagebox.showinfo("Run summary", "\n".join(summary_lines))
                self.log(f"Done. Wrote:\n  {target_out}\n  {not_found_out}")
        except Exception as e:
            self.log(f"ERROR: {e}")
            messagebox.showerror("Error", str(e))
        finally:
            self.run_btn_after_enable()

    def run_btn_after_enable(self):
        self.root.after(10, lambda: self.run_btn.config(state="normal"))

def main():
    root = tk.Tk()
    app = App(root)
    root.geometry("1200x760")
    root.mainloop()

if __name__ == "__main__":
    main()