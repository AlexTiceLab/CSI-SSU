#!/usr/bin/env python3
"""
Summarize taxonomic classification results from SQLite database.
"""

import sys
import os
import sqlite3
import csv
from collections import defaultdict, Counter

def generate_summary_csvs(db_path, taxonomy_summary_file, sequence_classifications_file, rank_counts_file, input_basename):
    """
    Generate CSV summary files from taxonomic classification database.
    
    Args:
        db_path: Path to SQLite results database
        taxonomy_summary_file: Path to output taxonomy summary CSV
        sequence_classifications_file: Path to output sequence classifications CSV
        rank_counts_file: Path to output rank counts CSV
        input_basename: Basename of input file for metadata
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(taxonomy_summary_file), exist_ok=True)
    
    try:
        # Connect to SQLite database
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Get basic statistics
        cursor.execute('SELECT COUNT(*) FROM placement_names')
        total_sequences = cursor.fetchone()[0]
        
        cursor.execute('SELECT COUNT(*) FROM taxa')
        total_taxa = cursor.fetchone()[0]
        
        # Get available taxonomic ranks from the database
        cursor.execute('SELECT DISTINCT want_rank FROM multiclass')
        available_ranks_set = set([row[0] for row in cursor.fetchall()])
        
        # Define taxonomic hierarchy order (from broad to specific)
        taxonomic_hierarchy = [
            'supergroup', 'division', 'subdivision', 'class', 'order', 'family', 'genus'
        ]
        
        # Use only the ranks that exist in the database, in proper hierarchical order
        ranks_of_interest = [rank for rank in taxonomic_hierarchy if rank in available_ranks_set]
        
        # Fall back to detected ranks if none match the hierarchy
        if not ranks_of_interest:
            ranks_of_interest = sorted(list(available_ranks_set))
        
        print(f"Detected taxonomic ranks: {ranks_of_interest}")
        
        # Get taxonomic classifications at different ranks
        classifications = {}
        
        for rank in ranks_of_interest:
            cursor.execute('''
            SELECT pn.name, t.tax_name, mc.likelihood
            FROM multiclass mc
            JOIN placement_names pn ON mc.placement_id = pn.placement_id
            JOIN taxa t ON mc.tax_id = t.tax_id
            WHERE mc.want_rank = ? AND mc.rank = ?
            ORDER BY pn.name, mc.likelihood DESC
            ''', (rank, rank))
            
            classifications[rank] = cursor.fetchall()
        
        # Get summary counts by taxonomic group at each rank
        taxonomic_counts = {}
        for rank in ranks_of_interest:
            cursor.execute('''
            SELECT t.tax_name, COUNT(DISTINCT pn.name) as seq_count, AVG(mc.likelihood) as avg_likelihood
            FROM multiclass mc
            JOIN placement_names pn ON mc.placement_id = pn.placement_id
            JOIN taxa t ON mc.tax_id = t.tax_id
            WHERE mc.want_rank = ? AND mc.rank = ?
            GROUP BY t.tax_name
            ORDER BY seq_count DESC, avg_likelihood DESC
            ''', (rank, rank))
            
            taxonomic_counts[rank] = cursor.fetchall()
        
        # Generate CSV files for downstream analysis
        
        # 1. Taxonomy Summary CSV - aggregated results by taxonomic rank and taxon
        with open(taxonomy_summary_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['rank', 'taxon', 'sequence_count', 'avg_likelihood'])
            
            for rank in ranks_of_interest:
                if rank in taxonomic_counts and taxonomic_counts[rank]:
                    for taxon, count, avg_likelihood in taxonomic_counts[rank]:
                        writer.writerow([rank, taxon, count, avg_likelihood])
        
        # 2. Sequence Classifications CSV - detailed lineage for each sequence
        with open(sequence_classifications_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Build header dynamically based on available ranks
            header = ['sequence_name']
            for rank in ranks_of_interest:
                header.extend([rank, f'{rank}_likelihood'])
            writer.writerow(header)
            
            # Get all sequence names
            cursor.execute('SELECT DISTINCT name FROM placement_names ORDER BY name')
            sequence_names = [row[0] for row in cursor.fetchall()]
            
            for seq_name in sequence_names:
                # Get full taxonomic lineage for this sequence
                lineage_dict = {}
                # Build placeholders for the IN clause
                placeholders = ','.join('?' * len(ranks_of_interest))
                
                # Build ORDER BY clause to match hierarchical order
                order_cases = []
                for i, rank in enumerate(ranks_of_interest):
                    order_cases.append(f"WHEN '{rank}' THEN {i}")
                order_clause = f"CASE mc.want_rank {' '.join(order_cases)} ELSE 999 END"
                
                # Modified query: Select only the highest likelihood classification for each rank
                # Uses subquery to get max likelihood per rank, then joins to get the corresponding taxon
                cursor.execute(f'''
                SELECT mc.want_rank, t.tax_name, mc.likelihood
                FROM multiclass mc
                JOIN placement_names pn ON mc.placement_id = pn.placement_id
                JOIN taxa t ON mc.tax_id = t.tax_id
                WHERE pn.name = ? AND mc.want_rank IN ({placeholders})
                AND mc.likelihood = (
                    SELECT MAX(mc2.likelihood)
                    FROM multiclass mc2
                    JOIN placement_names pn2 ON mc2.placement_id = pn2.placement_id
                    WHERE pn2.name = pn.name AND mc2.want_rank = mc.want_rank
                )
                ORDER BY {order_clause}
                ''', [seq_name] + ranks_of_interest)
                
                lineage = cursor.fetchall()
                # Keep only the highest likelihood for each rank (first occurrence if tied)
                for rank, taxon, likelihood in lineage:
                    if rank not in lineage_dict:
                        lineage_dict[rank] = (taxon, likelihood)
                
                # Write row with all taxonomic levels
                row = [seq_name]
                for rank in ranks_of_interest:
                    if rank in lineage_dict:
                        row.extend([lineage_dict[rank][0], lineage_dict[rank][1]])
                    else:
                        row.extend(['', ''])
                writer.writerow(row)
        
        # 3. Rank Counts CSV - summary statistics by rank
        with open(rank_counts_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['rank', 'unique_taxa', 'total_sequences'])
            
            for rank in ranks_of_interest:
                if rank in taxonomic_counts and taxonomic_counts[rank]:
                    unique_taxa = len(taxonomic_counts[rank])
                    classified_seqs = sum([count[1] for count in taxonomic_counts[rank]])
                    writer.writerow([rank, unique_taxa, classified_seqs])
        
        conn.close()
        
        print("CSV files generated successfully:")
        print("- Taxonomy summary: " + taxonomy_summary_file)
        print("- Sequence classifications: " + sequence_classifications_file) 
        print("- Rank counts: " + rank_counts_file)
        
    except Exception as e:
        # Write error to both log and output files for debugging
        error_msg = "Error processing results database: " + str(e)
        
        # Create empty CSV files with headers to prevent pipeline failure
        with open(taxonomy_summary_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['rank', 'taxon', 'sequence_count', 'avg_likelihood'])
        
        with open(sequence_classifications_file, 'w', newline='') as f:
            writer = csv.writer(f)
            # Use basic header if error occurs
            basic_header = ['sequence_name', 'rank', 'taxon', 'likelihood']
            writer.writerow(basic_header)
        
        with open(rank_counts_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['rank', 'unique_taxa', 'total_sequences'])
        
        raise Exception(error_msg)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python summarize_results.py <db_path> <taxonomy_summary_csv> <sequence_classifications_csv> <rank_counts_csv> <input_basename>")
        sys.exit(1)
    
    db_path = sys.argv[1]
    taxonomy_summary_file = sys.argv[2]
    sequence_classifications_file = sys.argv[3]
    rank_counts_file = sys.argv[4]
    input_basename = sys.argv[5]
    
    generate_summary_csvs(db_path, taxonomy_summary_file, sequence_classifications_file, rank_counts_file, input_basename)