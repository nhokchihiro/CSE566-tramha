import os
import glob
import csv
from cyvcf2 import VCF
import argparse

output_rows = []
control_rows = []

# Extract from vcf files only columns: sample, tool, chrom, start, end, svtype

def parse_vcf(path, sample, tool, output_rows):
    try:
        vcf = VCF(path)
        for var in vcf:
            chrom = var.CHROM
            start = var.start
            end = var.INFO.get('END')
            svtype = var.INFO.get('SVTYPE')
            if not end or not svtype:
                continue
            output_rows.append([sample, tool, chrom, start, end, svtype])
    except Exception as e:
        print(f"Error parsing {path}: {e}")

def parse_breakdancer(path, sample, output_rows):
    try:
        with open(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                chrom1 = parts[0]
                pos1 = int(parts[1])
                chrom2 = parts[3]
                pos2 = int(parts[4])
                svtype = parts[6]
                size = parts[7]
                output_rows.append([sample, 'breakdancer', chrom1, pos1, pos2, svtype])
    except Exception as e:
        print(f"Error parsing {path}: {e}")

def parse_lumpy(path, sample, output_rows):
    try:
        with open(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    print(f"Skipping malformed line: {line.strip()}")
                    continue
                chrom = parts[0]
                pos = int(parts[1])
                alt = parts[4]
                info = parts[7]
                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                    else:
                        info_dict[item] = True
                svtype = info_dict.get('SVTYPE', '')
                end = info_dict.get('END', pos)
                output_rows.append([sample, 'lumpy', chrom, pos, end, svtype])
    except Exception as e:
        print(f"Error parsing {path}: {e}")

def parse_manta(path, sample, output_rows):
    try:
        with open(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                end = None
                svtype = None
                for field in parts[7].split(';'):
                    if field.startswith('END='):
                        end = int(field.split('=')[1])
                    if field.startswith('SVTYPE='):
                        svtype = field.split('=')[1]
                output_rows.append([sample, 'manta', chrom, start, end, svtype])
    except Exception as e:
        print(f"Error parsing {path}: {e}")

def parse_delly(path, sample, output_rows):
    try:
        with open(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    print(f"Skipping malformed line: {line.strip()}")
                    continue
                chrom = parts[0]
                start = int(parts[1])
                info_fields = {}
                for field in parts[7].split(';'):
                    if '=' in field:
                        key, value = field.split('=', 1)
                        info_fields[key] = value
                    else:
                        info_fields[field] = True
                end = int(info_fields.get('END', start))
                svtype = info_fields.get('SVTYPE', '<DEL>')
                output_rows.append([sample, 'delly', chrom, start, end, svtype])
    except Exception as e:
        print(f"Error parsing {path}: {e}")

def parse_control_vcf(path, sample, control_rows):
    try:
        vcf = VCF(path)
        for var in vcf:
            chrom = var.CHROM
            start = var.start
            end = var.INFO.get('END')
            svtype = var.INFO.get('SVTYPE')
            if not end or not svtype:
                continue
            if svtype in ['INS', 'DEL']:
                control_rows.append([sample, 'SURVIVOR_sim', chrom, start, end, svtype])
    except Exception as e:
        print(f"Error parsing {path}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Extract SV calls from VCF files")
    parser.add_argument("--base_dir", required=True, help="Directory containing sample VCF files")
    parser.add_argument("--output_dir", default=".", help="Directory to write output files")
    args = parser.parse_args()

    for sample_folder in os.listdir(args.base_dir):
        sample_path = os.path.join(args.base_dir, sample_folder)
        if not os.path.isdir(sample_path):
            continue
        sample = sample_folder
        vcf_files = glob.glob(os.path.join(sample_path, f"{sample}.*.vcf"))

        for filepath in vcf_files:
            if 'simulated' in filepath.lower():
                parse_control_vcf(filepath, sample, control_rows)
                continue  # Skip further parsing for simulated files
            
            tool_name = os.path.basename(filepath).split('.')[1].lower()
            if tool_name == 'breakdancer':
                parse_breakdancer(filepath, sample, output_rows)
            elif tool_name == 'lumpy':
                parse_lumpy(filepath, sample, output_rows)
            elif tool_name == 'manta':
                parse_manta(filepath, sample, output_rows)
            elif tool_name == 'delly':
                parse_delly(filepath, sample, output_rows)
            else:
                parse_vcf(filepath, sample, tool_name, output_rows)

    with open(os.path.join(args.output_dir, "all_svcall.tsv"), "w", newline="") as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["sample", "tool", "chrom", "start", "end", "svtype"])
        writer.writerows(output_rows)

    with open(os.path.join(args.output_dir, "simulated_svcall.tsv"), "w", newline="") as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["sample", "tool", "chrom", "start", "end", "svtype"])
        writer.writerows(control_rows)

    print(f"Done! Extracted {len(output_rows)} SV entries to all_svcall.tsv and {len(control_rows)} control SV entries to simulated_svcall.tsv")


if __name__ == "__main__":
    main()
