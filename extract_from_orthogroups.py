#!/usr/bin/env python3
"""
CLI: extract_from_orthogroups.py

Reworked CLI to match GUI behavior with additional options:
 - Accepts a CSV queries file and a CSV column (or index)
 - Auto-detects whether the CSV column contains gene-prefixes and will enable
   prefix matching unless --no-auto-detect is passed and --match-gene-prefix is not set.
 - Optional --id-map (protein_id<TAB>gene_name) supported to annotate members if provided.
 - Supports target-species="All" (case-insensitive). If target == All, writes a single wide file:
     <out_prefix>_all_wide.tsv (wide format, one column per species EXCEPT source species)
 - If a specific target is chosen, file prefix includes the target species:
     <out_prefix>_<target>_target.tsv and <out_prefix>_<target>_not_found.tsv
 - Added --preview N to print first N unique IDs from the CSV column and exit.
 - Outputs no longer append any input CSV columns (GENE_ID, HD, JO, etc.)
"""
import argparse
import csv
import re
from collections import OrderedDict, defaultdict
import os
import sys

def normalize_protid_header(hdr):
    return hdr.split()[0] if hdr else ""

def gene_prefix_from_protid(protid):
    token = normalize_protid_header(protid)
    return re.split(r'[-|._:]', token)[0]

def sanitize_fname(s):
    return re.sub(r'\s+', '_', s)

def parse_orthogroups(path, log_fn=None):
    og_rows = OrderedDict()
    prot_to_og = {}
    prot_to_species = {}
    geneprefix_index = defaultdict(set)
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        if len(header) < 2:
            raise ValueError("Orthogroups.tsv header looks wrong.")
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

def format_ids(ids, id_map=None):
    if id_map:
        return "|".join([f"{i}|{id_map.get(i,'')}" for i in ids])
    return "|".join(ids)

def load_id_map(path):
    m = {}
    if not path:
        return m
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                m[parts[0]] = parts[1]
    return m

def preview_ids(rows_with_vals, n=50):
    seen = set()
    uniq = []
    for _, v in rows_with_vals:
        if not v:
            continue
        if v not in seen:
            seen.add(v)
            uniq.append(v)
        if len(uniq) >= n:
            break
    return uniq, len(seen)

def parse_keep_columns_arg(arg, header):
    """
    arg: comma-separated list, or "all"
    returns list of header names to keep (in order)
    """
    if not arg:
        return header[:]  # default all
    if arg.strip().lower() == "all":
        return header[:]
    parts = [p.strip() for p in arg.split(",") if p.strip()]
    keep = []
    for p in parts:
        if p in header:
            keep.append(p)
        else:
            # allow numeric indices
            try:
                idx = int(p)
                if 0 <= idx < len(header):
                    keep.append(header[idx])
                else:
                    print(f"Warning: keep-column index {p} out of range, ignored", file=sys.stderr)
            except Exception:
                print(f"Warning: keep-column '{p}' not in header and not an index, ignored", file=sys.stderr)
    return keep

def main():
    parser = argparse.ArgumentParser(description="Extract orthologs from OrthoFinder (CSV-support CLI)")
    parser.add_argument("--orthogroups", required=True)
    parser.add_argument("--csv", required=True, help="CSV file containing queries")
    parser.add_argument("--csv-column", required=True, help="CSV column name or 0-based index containing query IDs")
    parser.add_argument("--source-species", required=True, help="Species name in Orthogroups header (CSV IDs come from)")
    parser.add_argument("--target-species", required=True, help='Target species name in Orthogroups header to extract orthologs for OR "All"')
    parser.add_argument("--out-prefix", default="orthologs")
    parser.add_argument("--match-gene-prefix", action="store_true", help="Force match gene-prefix")
    parser.add_argument("--no-auto-detect", action="store_true", help="Disable auto-detection of CSV ID style")
    parser.add_argument("--id-map", help="optional TSV protein_id<TAB>gene_name to annotate output")
    parser.add_argument("--preview", type=int, default=0, help="If >0, print first N unique IDs and exit")
    # keep-columns option retained for compatibility but ignored (CLI will not append input CSV columns)
    parser.add_argument("--keep-columns", default="none", help='DEPRECATED: input CSV columns will not be appended to outputs')
    args = parser.parse_args()

    if not os.path.isfile(args.orthogroups):
        print("ERROR: orthogroups file not found:", args.orthogroups); sys.exit(1)
    if not os.path.isfile(args.csv):
        print("ERROR: csv file not found:", args.csv); sys.exit(1)

    species, og_rows, prot_to_og, prot_to_species, geneprefix_index = parse_orthogroups(args.orthogroups, log_fn=print)

    if args.source_species not in species:
        print("ERROR: source-species not one of orthogroups species. Available species:")
        print("\n".join(species)); sys.exit(1)

    target_is_all = (str(args.target_species).strip().lower() == "all")
    if not target_is_all and args.target_species not in species:
        print("ERROR: target-species not one of orthogroups species. Available species:")
        print("\n".join(species)); sys.exit(1)

    header, rows, delim = read_csv_rows(args.csv)
    if not header:
        print("ERROR: CSV empty or malformed"); sys.exit(1)

    rows_with_vals = extract_column_from_rows(rows, header, args.csv_column)

    # preview mode
    if args.preview and args.preview > 0:
        sample, total_unique = preview_ids(rows_with_vals, n=args.preview)
        print(f"Preview first {len(sample)} unique IDs (total unique in column: {total_unique}):")
        for s in sample:
            print(s)
        sys.exit(0)

    # NOTE: keep-columns option intentionally ignored — outputs will not include input CSV columns
    keep_cols = []  # explicit: do not append any input CSV columns to outputs

    # detect id style unless disabled or forced
    match_prefix = args.match_gene_prefix
    if not args.no_auto_detect and not match_prefix:
        detect = detect_column_id_style_from_rows(rows_with_vals, sample_n=200)
        if detect.get('mostly_prefix'):
            match_prefix = True
            print(f"Auto-detected CSV column '{args.csv_column}' as gene-prefix style (prefix_frac={detect['prefix_frac']:.2f}); enabling prefix matching.")
        else:
            print(f"CSV detection: protein_frac={detect['protein_frac']:.2f}, prefix_frac={detect['prefix_frac']:.2f} (no auto-change).")

    id_map = load_id_map(args.id_map) if args.id_map else {}

    total_queries = len(rows_with_vals)

    # All target: single wide file (exclude source species from wide columns)
    if target_is_all:
        out_fname = f"{sanitize_fname(args.out_prefix)}_all_wide.tsv"
        not_found_fname = f"{sanitize_fname(args.out_prefix)}_all_not_found.tsv"
        count_not_found = 0
        count_matched = 0
        count_matched_og = 0
        with open(out_fname, 'w', newline='') as wfh, open(not_found_fname, 'w', newline='') as nfh:
            writer = csv.writer(wfh, delimiter="\t")
            n_writer = csv.writer(nfh, delimiter="\t")
            species_excluding_source = [s for s in species if s != args.source_species]
            # header: no input CSV columns appended
            writer.writerow(["query_id", "orthogroup"] + species_excluding_source)
            n_writer.writerow(["query_id", "reason"])

            for i, (row, q) in enumerate(rows_with_vals, start=1):
                if i % 200 == 0:
                    print(f"Processing {i}/{total_queries}...")
                q = (q or "").strip()
                matched_proteins = set()
                if q and q in prot_to_og and prot_to_species.get(q) == args.source_species:
                    matched_proteins.add(q)
                normq = normalize_protid_header(q) if q else ""
                if normq and normq in prot_to_og and prot_to_species.get(normq) == args.source_species:
                    matched_proteins.add(normq)
                if match_prefix and q:
                    candidates = geneprefix_index.get(q, set())
                    for pid in candidates:
                        if prot_to_species.get(pid) == args.source_species:
                            matched_proteins.add(pid)
                    qtok = re.split(r'[-|._:]', normq)[0] if normq else ""
                    if qtok and qtok != q:
                        candidates2 = geneprefix_index.get(qtok, set())
                        for pid in candidates2:
                            if prot_to_species.get(pid) == args.source_species:
                                matched_proteins.add(pid)
                if not matched_proteins:
                    count_not_found += 1
                    n_writer.writerow([q, "NOT_FOUND"])
                    writer.writerow([q, ""] + [""]*len(species_excluding_source))
                    continue
                ogs = {prot_to_og[p] for p in matched_proteins if p in prot_to_og}
                if not ogs:
                    count_not_found += 1
                    n_writer.writerow([q, "NOT_FOUND"])
                    writer.writerow([q, ""] + [""]*len(species_excluding_source))
                    continue
                count_matched += 1
                count_matched_og += len(ogs)
                for og in sorted(ogs):
                    members_row = og_rows.get(og)
                    if members_row is None:
                        wide_cells = [""] * len(species_excluding_source)
                    else:
                        wide_cells = []
                        for sp in species_excluding_source:
                            try:
                                idx_sp = species.index(sp)
                                cell = members_row[idx_sp] if idx_sp < len(members_row) else ""
                            except ValueError:
                                cell = ""
                            ids = parse_members_cell(cell)
                            wide_cells.append(format_ids(ids, id_map=id_map))
                    writer.writerow([q, og] + wide_cells)
        summary_lines = [
            f"Total queries processed: {total_queries}",
            f"Queries with any OG match in source species: {count_matched}",
            f"Total OG-rows written (sum over matching OGs): {count_matched_og}",
            f"Queries NOT_FOUND: {count_not_found}"
        ]
        print("\n".join(summary_lines))
        print("Wrote:", out_fname, not_found_fname)
        sys.exit(0)

    # specific target: one file + not_found
    safe_target = sanitize_fname(args.target_species)
    prefix = f"{sanitize_fname(args.out_prefix)}_{safe_target}"
    target_out = f"{prefix}_target.tsv"
    not_found_out = f"{prefix}_not_found.tsv"

    count_not_found = 0
    count_matched = 0
    count_matched_target_nonempty = 0
    count_matched_og = 0

    with open(target_out, 'w', newline='') as tfh, open(not_found_out, 'w', newline='') as nfh:
        t_writer = csv.writer(tfh, delimiter="\t")
        n_writer = csv.writer(nfh, delimiter="\t")
        # header: no input CSV columns appended
        t_writer.writerow(["query_id", "orthogroup", f"{args.target_species}_members"])
        n_writer.writerow(["query_id", "reason"])
        for i, (row, q) in enumerate(rows_with_vals, start=1):
            if i % 200 == 0:
                print(f"Processing {i}/{total_queries}...")
            q = (q or "").strip()
            matched_proteins = set()
            if q and q in prot_to_og and prot_to_species.get(q) == args.source_species:
                matched_proteins.add(q)
            normq = normalize_protid_header(q) if q else ""
            if normq and normq in prot_to_og and prot_to_species.get(normq) == args.source_species:
                matched_proteins.add(normq)
            if match_prefix and q:
                candidates = geneprefix_index.get(q, set())
                for pid in candidates:
                    if prot_to_species.get(pid) == args.source_species:
                        matched_proteins.add(pid)
                qtok = re.split(r'[-|._:]', normq)[0] if normq else ""
                if qtok and qtok != q:
                    candidates2 = geneprefix_index.get(qtok, set())
                    for pid in candidates2:
                        if prot_to_species.get(pid) == args.source_species:
                            matched_proteins.add(pid)
            if not matched_proteins:
                count_not_found += 1
                n_writer.writerow([q, "NOT_FOUND"])
                t_writer.writerow([q, "", "NOT_FOUND"])
                continue
            ogs = {prot_to_og[p] for p in matched_proteins if p in prot_to_og}
            if not ogs:
                count_not_found += 1
                n_writer.writerow([q, "NOT_FOUND"])
                t_writer.writerow([q, "", "NOT_FOUND"])
                continue
            count_matched += 1
            count_matched_og += len(ogs)
            for og in sorted(ogs):
                members_row = og_rows.get(og)
                if members_row is None:
                    target_cell = ""
                else:
                    try:
                        idx_sp = species.index(args.target_species)
                        target_cell = members_row[idx_sp] if idx_sp < len(members_row) else ""
                    except ValueError:
                        target_cell = ""
                target_ids = parse_members_cell(target_cell)
                target_str = format_ids(target_ids, id_map=id_map)
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
    print("\n".join(summary_lines))
    print("Wrote:", target_out, not_found_out)

if __name__ == "__main__":
    main()
