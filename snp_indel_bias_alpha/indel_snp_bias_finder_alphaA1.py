import os
import sys
import gzip
import argparse
from pathlib import Path
from typing import List, Tuple, Set, Dict
import pysam

def get_vcf_file() -> str:
    """Get VCF file location from user."""
    if len(sys.argv) > 1:
        vcf_file = sys.argv[1]
    else:
        vcf_file = input("Enter the path to VCF file: ").strip()
    
    if not os.path.exists(vcf_file):
        raise FileNotFoundError(f"VCF file not found: {vcf_file}")
    
    return vcf_file

def create_or_find_index(vcf_file: str) -> str:
    """Find or create index file for VCF."""
    # Check if VCF is bgzipped
    if vcf_file.endswith('.gz'):
        index_file = vcf_file + '.tbi'
        
        # Check if index exists
        if not os.path.exists(index_file):
            print(f"Creating index file for {vcf_file}...")
            try:
                pysam.tabix_index(vcf_file, preset='vcf', force=True)
                print(f"Index created: {index_file}")
            except Exception as e:
                print(f"Warning: Could not create tabix index: {e}")
                print("Proceeding without index (slower processing)")
                return None
        else:
            print(f"Found existing index: {index_file}")
        
        return index_file
    else:
        print("VCF file is not bgzipped. For better performance, consider bgzipping it.")
        print("Proceeding without index...")
        return None

def is_indel(ref: str, alt: str) -> bool:
    """Check if variant is an indel based on REF and ALT."""
    # Split ALT by comma in case of multiple alleles
    alt_alleles = alt.split(',')
    
    # Check if REF or any ALT allele has length > 1
    if len(ref) > 1:
        return True
    
    for allele in alt_alleles:
        if len(allele) > 1:
            return True
    
    return False

def is_snp(ref: str, alt: str) -> bool:
    """Check if variant is a SNP (single nucleotide polymorphism)."""
    alt_alleles = alt.split(',')
    
    # SNP: REF is single base and all ALT alleles are single bases
    if len(ref) == 1:
        return all(len(allele) == 1 for allele in alt_alleles)
    
    return False

def find_snps_near_indels(vcf_file: str, window: int = 30) -> Tuple[List[str], List[Dict]]:
    """
    Find SNPs within specified window of indels.
    
    Args:
        vcf_file: Path to VCF file
        window: Base pair window around indels (default: 30)
    
    Returns:
        Tuple of (header_lines, snp_variants_near_indels)
    """
    header_lines = []
    all_variants = []
    indel_positions = {}  # {chrom: [positions]}
    
    # Determine if file is gzipped
    is_gzipped = vcf_file.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    print("Reading VCF file...")
    
    # First pass: Read all variants and identify indels
    with open_func(vcf_file, mode) as f:
        for line in f:
            if line.startswith('##'):
                header_lines.append(line.rstrip())
            elif line.startswith('#'):
                header_lines.append(line.rstrip())
                header_cols = line.rstrip().split('\t')
            else:
                fields = line.rstrip().split('\t')
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                
                variant_info = {
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'line': line.rstrip(),
                    'is_indel': is_indel(ref, alt),
                    'is_snp': is_snp(ref, alt)
                }
                
                all_variants.append(variant_info)
                
                # Store indel positions
                if variant_info['is_indel']:
                    if chrom not in indel_positions:
                        indel_positions[chrom] = []
                    indel_positions[chrom].append(pos)
    
    print(f"Total variants read: {len(all_variants)}")
    print(f"Total indels found: {sum(len(positions) for positions in indel_positions.values())}")
    
    # Second pass: Find SNPs near indels
    snps_near_indels = []
    snp_positions_used = set()  # To avoid duplicates
    
    print("Finding SNPs within 30 bases of indels...")
    
    for variant in all_variants:
        if variant['is_snp']:
            chrom = variant['chrom']
            pos = variant['pos']
            
            # Check if this chromosome has any indels
            if chrom in indel_positions:
                # Check if SNP is within window of any indel
                for indel_pos in indel_positions[chrom]:
                    if abs(pos - indel_pos) <= window and pos != indel_pos:
                        unique_key = (chrom, pos)
                        if unique_key not in snp_positions_used:
                            snps_near_indels.append(variant)
                            snp_positions_used.add(unique_key)
                        break  # No need to check other indels for this SNP
    
    print(f"SNPs found within {window} bases of indels: {len(snps_near_indels)}")
    
    return header_lines, snps_near_indels

def filter_clustered_variants(variants: List[Dict], window: int = 30, min_neighbors: int = 2) -> List[Dict]:
    """
    Filter variants to keep only those that have at least min_neighbors nearby variants
    within the specified window on the same chromosome.
    
    Args:
        variants: List of variant dictionaries
        window: Base pair window to check for neighbors (default: 30)
        min_neighbors: Minimum number of nearby variants required (default: 2)
    
    Returns:
        List of filtered variant dictionaries
    """
    print(f"\nApplying clustering filter (window={window}, min_neighbors={min_neighbors})...")
    
    # Group variants by chromosome
    variants_by_chrom = {}
    for variant in variants:
        chrom = variant['chrom']
        if chrom not in variants_by_chrom:
            variants_by_chrom[chrom] = []
        variants_by_chrom[chrom].append(variant)
    
    # Sort variants by position within each chromosome
    for chrom in variants_by_chrom:
        variants_by_chrom[chrom].sort(key=lambda x: x['pos'])
    
    # Find variants with sufficient neighbors
    filtered_variants = []
    
    for chrom, chrom_variants in variants_by_chrom.items():
        positions = [v['pos'] for v in chrom_variants]
        
        for i, variant in enumerate(chrom_variants):
            pos = variant['pos']
            neighbor_count = 0
            
            # Count neighbors within window
            # Check backwards
            j = i - 1
            while j >= 0 and positions[j] >= pos - window:
                if positions[j] != pos:  # Don't count itself
                    neighbor_count += 1
                j -= 1
            
            # Check forwards
            j = i + 1
            while j < len(positions) and positions[j] <= pos + window:
                if positions[j] != pos:  # Don't count itself
                    neighbor_count += 1
                j += 1
            
            # Keep variant if it has enough neighbors
            # Note: min_neighbors means we need at least that many OTHER variants nearby
            # So the total cluster size is min_neighbors + 1 (including the variant itself)
            if neighbor_count >= min_neighbors - 1:  # -1 because we're counting others, not including self
                filtered_variants.append(variant)
    
    print(f"Variants after clustering filter: {len(filtered_variants)}")
    print(f"Variants removed (isolated): {len(variants) - len(filtered_variants)}")
    
    return filtered_variants

def write_output_vcf(header_lines: List[str], variants: List[Dict], output_file: str):
    """Write filtered variants to output VCF file."""
    print(f"\nWriting output to: {output_file}")
    
    with open(output_file, 'w') as f:
        # Write header
        for header in header_lines:
            f.write(header + '\n')
        
        # Write variant lines
        for variant in variants:
            f.write(variant['line'] + '\n')
    
    print(f"Output VCF created successfully with {len(variants)} variants")

def main():
    parser = argparse.ArgumentParser(
        description='Find clustered SNPs within 30 bases of indels in VCF file'
    )
    parser.add_argument(
        'vcf_file',
        nargs='?',
        help='Path to input VCF file (optional, will prompt if not provided)'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output VCF file (default: input_filename_clustered_snps.vcf)',
        default=None
    )
    parser.add_argument(
        '-w', '--window',
        type=int,
        default=30,
        help='Base pair window around indels and for clustering (default: 30)'
    )
    parser.add_argument(
        '-m', '--min-neighbors',
        type=int,
        default=2,
        help='Minimum cluster size (including the variant itself, default: 2)'
    )
    
    args = parser.parse_args()
    
    try:
        # Get VCF file
        if args.vcf_file:
            vcf_file = args.vcf_file
            if not os.path.exists(vcf_file):
                raise FileNotFoundError(f"VCF file not found: {vcf_file}")
        else:
            vcf_file = get_vcf_file()
        
        print(f"Processing VCF file: {vcf_file}")
        print(f"Window size: {args.window} bases")
        print(f"Minimum cluster size: {args.min_neighbors} variants")
        
        # Create or find index
        create_or_find_index(vcf_file)
        
        # Step 1: Find SNPs near indels
        header_lines, snp_variants = find_snps_near_indels(vcf_file, window=args.window)
        
        # Step 2: Filter to keep only clustered variants
        # A variant is kept if it has at least (min_neighbors - 1) other variants within window
        clustered_variants = filter_clustered_variants(
            snp_variants, 
            window=args.window, 
            min_neighbors=args.min_neighbors
        )
        
        # Determine output filename
        if args.output:
            output_file = args.output
        else:
            base_name = os.path.splitext(os.path.basename(vcf_file))[0]
            if base_name.endswith('.vcf'):
                base_name = base_name[:-4]
            output_file = f"{base_name}_clustered_snps.vcf"
        
        # Write output
        write_output_vcf(header_lines, clustered_variants, output_file)
        
        print("\n" + "="*60)
        print("Processing complete!")
        print("="*60)
        print(f"Input file: {vcf_file}")
        print(f"Output file: {output_file}")
        print(f"Total SNPs near indels: {len(snp_variants)}")
        print(f"Clustered SNPs (final): {len(clustered_variants)}")
        print(f"Isolated SNPs removed: {len(snp_variants) - len(clustered_variants)}")
        print("="*60)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()