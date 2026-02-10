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

def categorize_indels_by_snp_proximity(vcf_file: str, window: int = 30) -> Tuple[List[str], Dict[str, List[Dict]], Dict[str, List[Dict]]]:
    """
    Categorize indels as biased (has nearby SNPs) or unbiased (no nearby SNPs).
    Also collect the SNPs near biased indels.
    
    Args:
        vcf_file: Path to VCF file
        window: Base pair window around indels (default: 30)
    
    Returns:
        Tuple of (header_lines, biased_data, unbiased_data)
        where biased_data and unbiased_data are dicts with 'indels' and 'snps' keys
    """
    header_lines = []
    all_variants = []
    
    # Determine if file is gzipped
    is_gzipped = vcf_file.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    print("Reading VCF file...")
    
    # Read all variants
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
    
    print(f"Total variants read: {len(all_variants)}")
    
    # Separate indels and SNPs
    indels = [v for v in all_variants if v['is_indel']]
    snps = [v for v in all_variants if v['is_snp']]
    
    print(f"Total indels found: {len(indels)}")
    print(f"Total SNPs found: {len(snps)}")
    
    # Group SNPs by chromosome for efficient lookup
    snps_by_chrom = {}
    for snp in snps:
        chrom = snp['chrom']
        if chrom not in snps_by_chrom:
            snps_by_chrom[chrom] = []
        snps_by_chrom[chrom].append(snp)
    
    # Sort SNPs by position within each chromosome
    for chrom in snps_by_chrom:
        snps_by_chrom[chrom].sort(key=lambda x: x['pos'])
    
    # Categorize indels and collect nearby SNPs
    biased_indels = []
    unbiased_indels = []
    biased_snps_set = set()  # Use set to avoid duplicates
    biased_snps = []
    
    print(f"\nCategorizing indels based on nearby SNPs (window={window})...")
    
    for indel in indels:
        chrom = indel['chrom']
        pos = indel['pos']
        has_nearby_snp = False
        
        # Check if this chromosome has SNPs
        if chrom in snps_by_chrom:
            # Find SNPs within window
            for snp in snps_by_chrom[chrom]:
                snp_pos = snp['pos']
                
                # Check if SNP is within window
                if abs(snp_pos - pos) <= window and snp_pos != pos:
                    has_nearby_snp = True
                    
                    # Add SNP to biased SNPs collection (avoid duplicates)
                    snp_key = (snp['chrom'], snp['pos'])
                    if snp_key not in biased_snps_set:
                        biased_snps_set.add(snp_key)
                        biased_snps.append(snp)
        
        # Categorize indel
        if has_nearby_snp:
            biased_indels.append(indel)
        else:
            unbiased_indels.append(indel)
    
    print(f"Biased indels (with nearby SNPs): {len(biased_indels)}")
    print(f"Unbiased indels (no nearby SNPs): {len(unbiased_indels)}")
    print(f"SNPs near biased indels: {len(biased_snps)}")
    
    # Prepare return data
    biased_data = {
        'indels': biased_indels,
        'snps': biased_snps
    }
    
    unbiased_data = {
        'indels': unbiased_indels,
        'snps': []  # No SNPs for unbiased dataset
    }
    
    return header_lines, biased_data, unbiased_data

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
    if len(variants) == 0:
        return []
    
    print(f"Applying clustering filter (window={window}, min_neighbors={min_neighbors})...")
    
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
            if neighbor_count >= min_neighbors - 1:
                filtered_variants.append(variant)
    
    print(f"Variants after clustering filter: {len(filtered_variants)}")
    print(f"Variants removed (isolated): {len(variants) - len(filtered_variants)}")
    
    return filtered_variants

def write_output_vcf(header_lines: List[str], variants: List[Dict], output_file: str):
    """Write filtered variants to output VCF file."""
    print(f"Writing output to: {output_file}")
    
    with open(output_file, 'w') as f:
        # Write header
        for header in header_lines:
            f.write(header + '\n')
        
        # Write variant lines (sorted by chromosome and position)
        # Sort for better output organization
        sorted_variants = sorted(variants, key=lambda x: (x['chrom'], x['pos']))
        
        for variant in sorted_variants:
            f.write(variant['line'] + '\n')
    
    print(f"Output VCF created successfully with {len(variants)} variants")

def main():
    parser = argparse.ArgumentParser(
        description='Segregate indels into biased (with nearby SNPs) and unbiased (without nearby SNPs) datasets'
    )
    parser.add_argument(
        'vcf_file',
        nargs='?',
        help='Path to input VCF file (optional, will prompt if not provided)'
    )
    parser.add_argument(
        '-p', '--prefix',
        help='Output file prefix (default: input filename without extension)',
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
        help='Minimum cluster size for SNPs in biased dataset (default: 2)'
    )
    parser.add_argument(
        '--no-clustering',
        action='store_true',
        help='Skip clustering filter for biased SNPs'
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
        if not args.no_clustering:
            print(f"Minimum cluster size for biased SNPs: {args.min_neighbors} variants")
        
        # Create or find index
        create_or_find_index(vcf_file)
        
        # Step 1: Categorize indels by SNP proximity
        header_lines, biased_data, unbiased_data = categorize_indels_by_snp_proximity(
            vcf_file, 
            window=args.window
        )
        
        # Step 2: Apply clustering filter to biased SNPs if requested
        if not args.no_clustering and len(biased_data['snps']) > 0:
            print("\nFiltering biased SNPs for clustering...")
            clustered_snps = filter_clustered_variants(
                biased_data['snps'],
                window=args.window,
                min_neighbors=args.min_neighbors
            )
            biased_data['snps'] = clustered_snps
        
        # Determine output prefix
        if args.prefix:
            prefix = args.prefix
        else:
            base_name = os.path.splitext(os.path.basename(vcf_file))[0]
            if base_name.endswith('.vcf'):
                base_name = base_name[:-4]
            prefix = base_name
        
        # Create output filenames
        biased_output = f"{prefix}_biased.vcf"
        unbiased_output = f"{prefix}_unbiased.vcf"
        
        # Combine indels and SNPs for biased dataset
        biased_variants = biased_data['indels'] + biased_data['snps']
        unbiased_variants = unbiased_data['indels']
        
        # Write outputs
        print("\n" + "="*60)
        print("Writing output files...")
        print("="*60)
        
        write_output_vcf(header_lines, biased_variants, biased_output)
        print()
        write_output_vcf(header_lines, unbiased_variants, unbiased_output)
        
        # Final summary
        print("\n" + "="*60)
        print("Processing complete!")
        print("="*60)
        print(f"Input file: {vcf_file}")
        print(f"Output prefix: {prefix}")
        print()
        print("BIASED DATASET (indels with nearby SNPs):")
        print(f"  File: {biased_output}")
        print(f"  Indels: {len(biased_data['indels'])}")
        print(f"  SNPs: {len(biased_data['snps'])}")
        print(f"  Total variants: {len(biased_variants)}")
        print()
        print("UNBIASED DATASET (indels without nearby SNPs):")
        print(f"  File: {unbiased_output}")
        print(f"  Indels: {len(unbiased_data['indels'])}")
        print(f"  SNPs: {len(unbiased_data['snps'])}")
        print(f"  Total variants: {len(unbiased_variants)}")
        print("="*60)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()