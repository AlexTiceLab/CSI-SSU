#!/usr/bin/env python3
"""
Plot phylogenetic tree with placed sequences highlighted using ete3.
"""

import sys
import os
import sqlite3

def get_sequence_lengths(fasta_file):
    """
    Read sequence lengths from a FASTA file.
    
    Args:
        fasta_file: Path to FASTA file
    
    Returns:
        dict: Dictionary mapping sequence names to their lengths
    """
    sequence_lengths = {}
    
    if not fasta_file or not os.path.exists(fasta_file):
        print(f"Warning: FASTA file not found: {fasta_file}")
        return sequence_lengths
    
    try:
        current_name = None
        current_length = 0
        
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence length
                    if current_name is not None:
                        sequence_lengths[current_name] = current_length
                    
                    # Start new sequence
                    current_name = line[1:]  # Remove '>' character
                    current_length = 0
                elif line and current_name is not None:
                    # Count sequence characters (ignore gaps and ambiguous characters)
                    current_length += len([c for c in line if c.isalpha() and c not in '-Nn'])
            
            # Don't forget the last sequence
            if current_name is not None:
                sequence_lengths[current_name] = current_length
        
        print(f"Read lengths for {len(sequence_lengths)} sequences from {fasta_file}")
        return sequence_lengths
        
    except Exception as e:
        print(f"Warning: Could not read sequence lengths from {fasta_file}: {e}")
        return sequence_lengths

def root_tree_by_supergroup(tree, db_path, supergroup_of_interest):
    """
    Root the tree between the supergroup of interest and all other taxa.
    
    Args:
        tree: ete3 Tree object
        db_path: Path to SQLite results database
        supergroup_of_interest: Name of supergroup to root by
    
    Returns:
        Tree: Rooted tree object
    """
    if not db_path or not os.path.exists(db_path):
        print("Warning: No database provided. Cannot root tree by supergroup.")
        return tree
    
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Get all reference sequences and their taxonomic classifications
        cursor.execute('''
        SELECT DISTINCT t.tax_name, t.rank
        FROM taxa t
        WHERE t.rank = 'supergroup'
        ''')
        
        supergroups = cursor.fetchall()
        print(f"Available supergroups in database: {[sg[0] for sg in supergroups]}")
        
        # Get sequences that belong to the supergroup of interest
        cursor.execute('''
        SELECT DISTINCT pn.name
        FROM placement_names pn
        JOIN multiclass mc ON pn.placement_id = mc.placement_id
        JOIN taxa t ON mc.tax_id = t.tax_id
        WHERE t.tax_name = ? AND t.rank = 'supergroup'
        ''', (supergroup_of_interest,))
        
        supergroup_sequences = {row[0] for row in cursor.fetchall()}
        print(f"Found {len(supergroup_sequences)} sequences from {supergroup_of_interest}")
        
        conn.close()
        
        # Get all leaf names from tree
        all_leaves = {leaf.name for leaf in tree.get_leaves() if leaf.name}
        
        # Find leaves NOT in the supergroup of interest
        non_supergroup_leaves = all_leaves - supergroup_sequences
        print(f"Found {len(non_supergroup_leaves)} sequences NOT from {supergroup_of_interest}")
        
        if len(supergroup_sequences) == 0:
            print(f"Warning: No sequences found for supergroup '{supergroup_of_interest}'. Cannot root tree.")
            return tree
        
        if len(non_supergroup_leaves) == 0:
            print(f"Warning: All sequences belong to {supergroup_of_interest}. Cannot root tree.")
            return tree
        
        # Find the most recent common ancestor (MRCA) of non-supergroup sequences
        non_supergroup_nodes = []
        for leaf in tree.get_leaves():
            if leaf.name and leaf.name in non_supergroup_leaves:
                non_supergroup_nodes.append(leaf)
        
        if len(non_supergroup_nodes) < 2:
            print("Warning: Need at least 2 non-supergroup sequences to root tree.")
            return tree
        
        # Get MRCA of non-supergroup sequences to use as outgroup
        outgroup_mrca = tree.get_common_ancestor(non_supergroup_nodes)
        
        # Root the tree using this outgroup
        tree.set_outgroup(outgroup_mrca)
        
        tree.ladderize(direction=0) 
        
        print(f"Successfully rooted tree with {supergroup_of_interest} as ingroup and {len(non_supergroup_leaves)} sequences as outgroup")
        
        return tree
        
    except Exception as e:
        print(f"Warning: Could not root tree by supergroup: {e}")
        print("Continuing with unrooted tree.")
        return tree

def get_taxonomic_hierarchy(db_path, reference_package_dir=None):
    """
    Get taxonomic hierarchy for all sequences.
    Placed sequences: [supergroup|division|order|family]
    Reference sequences: [supergroup|division|order|species]
    Uses SQLite database for placed sequences and reference package for reference sequences.
    
    Args:
        db_path: Path to SQLite results database
        reference_package_dir: Path to reference package directory containing {supergroup}.tax and {supergroup}.seq_info
    
    Returns:
        dict: Dictionary mapping sequence names to their taxonomic hierarchy string
    """
    taxonomic_info = {}
    placed_tax_info = {}  # Initialize here to avoid scoping issues
    
    # Get taxonomy for placed sequences from SQLite database
    if db_path and os.path.exists(db_path):
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Get all sequences and their complete taxonomic hierarchy
            cursor.execute('''
            SELECT DISTINCT pn.name, t.tax_name, t.rank
            FROM placement_names pn
            JOIN multiclass mc ON pn.placement_id = mc.placement_id
            JOIN taxa t ON mc.tax_id = t.tax_id
            WHERE t.rank IN ('supergroup', 'division', 'order', 'family', 'genus')
            ORDER BY pn.name, 
                     CASE t.rank 
                         WHEN 'supergroup' THEN 1
                         WHEN 'division' THEN 2
                         WHEN 'order' THEN 3
                         WHEN 'family' THEN 4
                         WHEN 'genus' THEN 5
                     END
            ''')
            
            rows = cursor.fetchall()
            
            # Build hierarchy for placed sequences
            for name, taxon, rank in rows:
                if name not in placed_tax_info:
                    placed_tax_info[name] = {'supergroup': '', 'division': '', 'order': '', 'family': '', 'genus': ''}
                placed_tax_info[name][rank] = taxon
            
            # Convert to hierarchy strings for placed sequences (supergroup|division|order|family)
            for name, taxa in placed_tax_info.items():
                hierarchy = f"{taxa['supergroup']}|{taxa['division']}|{taxa['order']}|{taxa['family']}"
                # Clean up empty fields
                hierarchy = '|'.join([part if part else 'Uncertain' for part in hierarchy.split('|')])
                taxonomic_info[name] = hierarchy
            
            conn.close()
            print(f"Retrieved taxonomic hierarchy for {len(taxonomic_info)} placed sequences from database")
            
        except Exception as e:
            print(f"Warning: Could not read taxonomic hierarchy from database: {e}")
    
    # Get taxonomy for reference sequences from reference package
    if reference_package_dir and os.path.exists(reference_package_dir):
        # Get supergroup name from reference package directory
        supergroup_name = os.path.basename(reference_package_dir).replace('.refpkg', '')
        
        # Use unified naming convention: {supergroup}.tax
        taxonomy_file = os.path.join(reference_package_dir, f'{supergroup_name}.tax')
        
        if os.path.exists(taxonomy_file):
            try:
                import csv
                
                print(f"Reading taxonomy file: {taxonomy_file}")
                with open(taxonomy_file, 'r') as f:
                    reader = csv.DictReader(f)
                    
                    # Build mappings from tax_id to both tax_name, rank, and parent relationships
                    tax_id_to_name = {}
                    tax_id_to_rank = {}
                    tax_id_to_parent = {}
                    tax_id_to_hierarchy = {}
                    
                    for row in reader:
                        tax_id = row['tax_id']
                        tax_name = row['tax_name']
                        rank = row['rank']
                        parent_id = row['parent_id']
                        
                        # Store the name, rank, and parent for this tax_id
                        tax_id_to_name[tax_id] = tax_name
                        tax_id_to_rank[tax_id] = rank
                        tax_id_to_parent[tax_id] = parent_id
                        
                        # Build hierarchical path by following parent relationships
                        if tax_id not in tax_id_to_hierarchy:
                            tax_id_to_hierarchy[tax_id] = {
                                'supergroup': '',
                                'division': '', 
                                'order': '',
                                'genus': '',
                                'species': ''
                            }
                        
                        # Store the taxonomic names at each level
                        supergroup_id = row.get('supergroup', '')
                        division_id = row.get('division', '')
                        order_id = row.get('order', '')
                        genus_id = row.get('genus', '')
                        species_id = row.get('species', '')
                        
                        if supergroup_id:
                            tax_id_to_hierarchy[tax_id]['supergroup'] = supergroup_id
                        if division_id:
                            tax_id_to_hierarchy[tax_id]['division'] = division_id
                        if order_id:
                            tax_id_to_hierarchy[tax_id]['order'] = order_id
                        if genus_id:
                            tax_id_to_hierarchy[tax_id]['genus'] = genus_id
                        if species_id:
                            tax_id_to_hierarchy[tax_id]['species'] = species_id
                
                # Second pass: convert tax_ids to tax_names in the hierarchy
                # For isolates, we want to show the parent species name instead
                tax_id_to_info = {}
                for tax_id, hierarchy in tax_id_to_hierarchy.items():
                    # Special handling for species level - if this is an isolate, use parent species
                    species_name = ''
                    if hierarchy['species']:
                        species_tax_id = hierarchy['species']
                        if tax_id_to_rank.get(species_tax_id) == 'species':
                            species_name = tax_id_to_name.get(species_tax_id, '')
                        elif tax_id_to_rank.get(tax_id) == 'isolate' and tax_id_to_parent.get(tax_id):
                            # For isolates, get the parent species name
                            parent_id = tax_id_to_parent.get(tax_id)
                            if tax_id_to_rank.get(parent_id) == 'species':
                                species_name = tax_id_to_name.get(parent_id, '')
                    
                    tax_id_to_info[tax_id] = {
                        'supergroup': tax_id_to_name.get(hierarchy['supergroup'], '') if hierarchy['supergroup'] else '',
                        'division': tax_id_to_name.get(hierarchy['division'], '') if hierarchy['division'] else '',
                        'order': tax_id_to_name.get(hierarchy['order'], '') if hierarchy['order'] else '',
                        'genus': tax_id_to_name.get(hierarchy['genus'], '') if hierarchy['genus'] else '',
                        'species': species_name
                    }
                
                # Now read the seq_info file to map sequence names to tax_ids
                # Use unified naming convention: {supergroup}.seq_info
                seq_info_file = os.path.join(reference_package_dir, f'{supergroup_name}.seq_info')
                
                if os.path.exists(seq_info_file):
                    reference_count = 0
                    print(f"Reading seq_info file: {seq_info_file}")
                    with open(seq_info_file, 'r') as f:
                        reader = csv.DictReader(f)
                        
                        for row in reader:
                            seq_name = row.get('seqname', '')
                            tax_id = row.get('tax_id', '')
                            
                            # Skip if this sequence already has taxonomy from placed sequences
                            if seq_name in taxonomic_info:
                                continue
                                
                            if tax_id in tax_id_to_info:
                                taxa = tax_id_to_info[tax_id]
                                # For reference sequences: supergroup|division|order|species
                                hierarchy = f"{taxa['supergroup']}|{taxa['division']}|{taxa['order']}|{taxa['species']}"
                                # Clean up empty fields
                                hierarchy = '|'.join([part if part else 'Uncertain' for part in hierarchy.split('|')])
                                taxonomic_info[seq_name] = hierarchy
                                reference_count += 1
                    
                    print(f"Retrieved taxonomic hierarchy for {reference_count} reference sequences from reference package")
                    
            except Exception as e:
                print(f"Warning: Could not read taxonomic hierarchy from reference package: {e}")
    
    print(f"Total sequences with taxonomic hierarchy: {len(taxonomic_info)}")
    return taxonomic_info

def get_placed_sequences(db_path):
    """
    Get list of placed sequence names from the SQLite database with supergroup, division, order and family-level likelihoods.
    
    Args:
        db_path: Path to SQLite results database
    
    Returns:
        tuple: (placed_sequences dict, supergroup_likelihoods dict, division_likelihoods dict, order_likelihoods dict, family_likelihoods dict)
    """
    placed_sequences = {}
    supergroup_likelihoods = {}
    division_likelihoods = {}
    order_likelihoods = {}
    family_likelihoods = {}
    
    if os.path.exists(db_path):
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Get placed sequences with their best classification at genus/species level
            cursor.execute('''
            SELECT pn.name, t.tax_name, t.rank, mc.likelihood
            FROM placement_names pn
            JOIN multiclass mc ON pn.placement_id = mc.placement_id
            JOIN taxa t ON mc.tax_id = t.tax_id
            WHERE mc.want_rank IN ('genus', 'species') AND mc.rank IN ('genus', 'species')
            ORDER BY pn.name, mc.likelihood DESC
            ''')
            
            # Get supergroup-level likelihoods for placed sequences
            cursor.execute('''
            SELECT pn.name, mc.likelihood
            FROM placement_names pn
            JOIN multiclass mc ON pn.placement_id = mc.placement_id
            JOIN taxa t ON mc.tax_id = t.tax_id
            WHERE mc.want_rank = 'supergroup' AND mc.rank = 'supergroup'
            ORDER BY pn.name, mc.likelihood DESC
            ''')
            
            supergroup_rows = cursor.fetchall()
            for name, likelihood in supergroup_rows:
                if name not in supergroup_likelihoods:
                    supergroup_likelihoods[name] = likelihood
            
            # Get division-level likelihoods for placed sequences
            cursor.execute('''
            SELECT pn.name, mc.likelihood
            FROM placement_names pn
            JOIN multiclass mc ON pn.placement_id = mc.placement_id
            JOIN taxa t ON mc.tax_id = t.tax_id
            WHERE mc.want_rank = 'division' AND mc.rank = 'division'
            ORDER BY pn.name, mc.likelihood DESC
            ''')
            
            division_rows = cursor.fetchall()
            for name, likelihood in division_rows:
                if name not in division_likelihoods:
                    division_likelihoods[name] = likelihood
            
            # Also get order-level likelihoods for placed sequences
            cursor.execute('''
            SELECT pn.name, mc.likelihood
            FROM placement_names pn
            JOIN multiclass mc ON pn.placement_id = mc.placement_id
            JOIN taxa t ON mc.tax_id = t.tax_id
            WHERE mc.want_rank = 'order' AND mc.rank = 'order'
            ORDER BY pn.name, mc.likelihood DESC
            ''')
            
            order_rows = cursor.fetchall()
            for name, likelihood in order_rows:
                if name not in order_likelihoods:
                    order_likelihoods[name] = likelihood
            
            # Also get family-level likelihoods for placed sequences
            cursor.execute('''
            SELECT pn.name, mc.likelihood
            FROM placement_names pn
            JOIN multiclass mc ON pn.placement_id = mc.placement_id
            JOIN taxa t ON mc.tax_id = t.tax_id
            WHERE mc.want_rank = 'family' AND mc.rank = 'family'
            ORDER BY pn.name, mc.likelihood DESC
            ''')
            
            family_rows = cursor.fetchall()
            for name, likelihood in family_rows:
                if name not in family_likelihoods:
                    family_likelihoods[name] = likelihood
            
            rows = cursor.fetchall()
            
            for name, taxon, rank, likelihood in rows:
                if name not in placed_sequences:
                    placed_sequences[name] = f"{taxon} ({rank})"
            
            # Also get sequences that may not have genus/species level classification
            cursor.execute('SELECT DISTINCT name FROM placement_names')
            all_placed = cursor.fetchall()
            
            for (name,) in all_placed:
                if name not in placed_sequences:
                    placed_sequences[name] = "Uncertain classification"
            
            conn.close()
        except Exception as e:
            print(f"Warning: Could not read placement database: {e}")
            print("Will plot tree without taxonomic information.")
    
    return placed_sequences, supergroup_likelihoods, division_likelihoods, order_likelihoods, family_likelihoods

def plot_tree_ete3(tree_file, output_file, db_path=None, supergroup_of_interest=None, reference_package_dir=None, fasta_file=None):
    """
    Plot phylogenetic tree using ete3 with placed sequences in bold.
    
    Args:
        tree_file: Path to input tree file (Newick format)
        output_file: Path to output image file
        db_path: Optional path to SQLite database for getting placed sequences
        supergroup_of_interest: Optional supergroup to use for rooting the tree
        reference_package_dir: Optional path to reference package directory
        fasta_file: Optional path to FASTA file containing sequence lengths for placed sequences
    """
    import os
    # Set environment variables for headless rendering
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    os.environ['DISPLAY'] = ''
    
    from ete3 import Tree, TreeStyle, NodeStyle, TextFace, AttrFace
    
    print(f"Loading tree from: {tree_file}")
    
    # Load the tree - handle square bracket support values by converting them to proper Newick format
    try:
        # First attempt: try loading directly (works if no square brackets)
        tree = Tree(tree_file)
    except Exception as e:
        print(f"Initial tree loading failed, trying to convert support value format...")
        # Handle square bracket support values by converting to proper Newick format
        import tempfile
        import re
        
        with open(tree_file, 'r') as f:
            tree_content = f.read()
        
        # Convert square bracket support values [XX] to proper Newick format
        # In Newick format, support values should come after the closing parenthesis, before the colon
        # Pattern: ):0.123[95] should become )95:0.123
        def convert_support_values(match):
            branch_length = match.group(1)  # The 0.123 part (without colon)
            support_value = match.group(2)  # The 95 part
            return f"){support_value}:{branch_length}"
        
        # Apply conversion: ):branch_length[support] -> )support:branch_length
        # Match pattern like ):0.0077643[0] or ):1e-06[100]
        cleaned_content = re.sub(r'\):([0-9\.e\-\+]+)\[([0-9\.]+)\]', convert_support_values, tree_content)
        
        # Also handle cases where there's no branch length: )[support] -> )support
        cleaned_content = re.sub(r'\)\[([0-9\.]+)\]', r')\1', cleaned_content)
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tre', delete=False) as temp_file:
            temp_file.write(cleaned_content)
            temp_tree_file = temp_file.name
        
        try:
            tree = Tree(temp_tree_file)
            print("Successfully converted square bracket support values to proper Newick format")
            # Clean up temporary file
            import os
            os.unlink(temp_tree_file)
        except Exception as e2:
            print(f"Support value conversion failed, trying without support values: {e2}")
            # Fallback: remove support values entirely but keep branch lengths
            cleaned_content = re.sub(r'\[([0-9\.]+)\]', '', tree_content)
            with tempfile.NamedTemporaryFile(mode='w', suffix='.tre', delete=False) as temp_file2:
                temp_file2.write(cleaned_content)
                temp_tree_file2 = temp_file2.name
            
            try:
                tree = Tree(temp_tree_file2)
                print("Warning: Loaded tree without support values - all support values will show as 1.0")
                # Clean up temporary files
                import os
                os.unlink(temp_tree_file)
                os.unlink(temp_tree_file2)
            except Exception as e3:
                print(f"All parsing attempts failed: {e3}")
                # Clean up temporary files
                import os
                try:
                    os.unlink(temp_tree_file)
                    os.unlink(temp_tree_file2)
                except:
                    pass
                raise ImportError(f"ete3 cannot parse this tree format. Original error: {e}")
    
    # Root the tree by supergroup if specified
    if supergroup_of_interest:
        tree = root_tree_by_supergroup(tree, db_path, supergroup_of_interest)
    
    # Get placed sequences, taxonomic hierarchy, sequence lengths, and likelihoods
    placed_sequences = {}
    supergroup_likelihoods = {}
    division_likelihoods = {}
    order_likelihoods = {}
    family_likelihoods = {}
    taxonomic_hierarchy = {}
    sequence_lengths = {}
    
    # Read sequence lengths from FASTA file
    if fasta_file:
        sequence_lengths = get_sequence_lengths(fasta_file)
    
    if db_path and os.path.exists(db_path):
        placed_sequences, supergroup_likelihoods, division_likelihoods, order_likelihoods, family_likelihoods = get_placed_sequences(db_path)
        taxonomic_hierarchy = get_taxonomic_hierarchy(db_path, reference_package_dir)
        print(f"Found {len(placed_sequences)} placed sequences with taxonomic information")
    else:
        print("No database provided. Looking for sequences with typical placement patterns.")
        # Still try to get reference taxonomy if available
        taxonomic_hierarchy = get_taxonomic_hierarchy(None, reference_package_dir)
        # Fallback: identify placed sequences by common patterns
        for leaf in tree:
            if leaf.name and ("|scaffold_" in leaf.name or "|contig_" in leaf.name or 
                             "P10K-" in leaf.name or "_out" in leaf.name):
                placed_sequences[leaf.name] = "Placed sequence (pattern match)"
    
    # Calculate dynamic figure dimensions - height scales with taxa
    num_taxa = len(list(tree.get_leaves()))
    
    # Height should scale with number of taxa for readability
    figure_height = max(600, num_taxa * 20)  # 20 pixels per taxon for good spacing
    # Width should be proportional but not too wide to avoid stretching
    figure_width = min(1200, max(800, figure_height * 0.8))  # Reasonable aspect ratio
    
    print(f"Creating ete3 figure for {num_taxa} taxa: {figure_width}x{figure_height} pixels")
    
    # Create tree style optimized for tall trees
    ts = TreeStyle()
    ts.show_leaf_name = False  # Disable built-in names since we're adding custom faces
    ts.show_branch_length = False  # Remove branch lengths from display
    ts.show_branch_support = False  # Disable automatic support display, we'll add custom ones
    ts.mode = "r"
    ts.tree_width = figure_width * 0.5  # Leave more room for labels
    ts.allow_face_overlap = False  # Prevent label overlap
    ts.branch_vertical_margin = 5  # Consistent spacing between branches
    ts.draw_aligned_faces_as_table = False
    ts.force_topology = False
    ts.margin_top = 10
    ts.margin_bottom = 10
    ts.margin_left = 10
    ts.margin_right = 10
    
    # Style for placed sequences
    placed_style = NodeStyle()
    placed_style["fgcolor"] = "black"
    placed_style["size"] = 0
    placed_style["shape"] = "circle"
    
    # Style for reference sequences
    reference_style = NodeStyle()
    reference_style["fgcolor"] = "black"
    reference_style["size"] = 0
    reference_style["shape"] = "circle"
    
    # Count placed sequences and set node styles
    formatted_count = 0
    
    # Style for internal nodes (remove blue circles)
    internal_style = NodeStyle()
    internal_style["size"] = 0
    
    # Apply styles to all nodes
    for node in tree.traverse():
        if node.is_leaf():
            if node.name:
                # Create display name with taxonomic hierarchy if available
                if node.name in taxonomic_hierarchy:
                    display_name = f"{node.name} [{taxonomic_hierarchy[node.name]}]"
                else:
                    display_name = node.name
                
                if node.name in placed_sequences:
                    # Add sequence length and order likelihood for placed sequences only
                    if node.name in sequence_lengths:
                        display_name = f"{node.name} ({sequence_lengths[node.name]}bp) [{taxonomic_hierarchy[node.name]}]"
                    
                    # Add supergroup, division, order and family likelihood if available
                    likelihood_parts = []
                    if node.name in supergroup_likelihoods:
                        likelihood_parts.append(f"S: {supergroup_likelihoods[node.name]:.3f}")
                    if node.name in division_likelihoods:
                        likelihood_parts.append(f"D: {division_likelihoods[node.name]:.3f}")
                    if node.name in order_likelihoods:
                        likelihood_parts.append(f"O: {order_likelihoods[node.name]:.3f}")
                    if node.name in family_likelihoods:
                        likelihood_parts.append(f"F: {family_likelihoods[node.name]:.3f}")
                    if likelihood_parts:
                        likelihood_str = f" ({', '.join(likelihood_parts)})"
                        display_name += likelihood_str
                    
                    # Create new style for placed sequences without support values
                    placed_no_support_style = NodeStyle()
                    placed_no_support_style["fgcolor"] = "black"
                    placed_no_support_style["size"] = 0
                    placed_no_support_style["shape"] = "circle"
                    node.set_style(placed_no_support_style)
                    # Add bold text face for placed sequences
                    name_face = TextFace(display_name, fgcolor="black", fsize=9, bold=True)
                    name_face.margin_left = 5
                    node.add_face(name_face, column=0, position="branch-right")
                    formatted_count += 1
                else:
                    node.set_style(reference_style)
                    # Add normal text face for reference sequences (no sequence length)
                    name_face = TextFace(display_name, fgcolor="black", fsize=8, bold=False)
                    name_face.margin_left = 5
                    node.add_face(name_face, column=0, position="branch-right")
        else:
            # Check if this internal node has any placed sequences as direct children
            has_placed_child = any(child.is_leaf() and child.name in placed_sequences for child in node.children)
            
            # Apply no-circle style to all internal nodes
            node.set_style(internal_style)
            
            # Add custom support value display only for nodes without placed children
            if not has_placed_child and hasattr(node, 'support') and node.support is not None and node.support >= 0:
                # Add support value as a custom text face for reference nodes only
                support_face = TextFace(str(int(node.support)), fgcolor="black", fsize=8)
                support_face.margin_left = 2
                support_face.margin_right = 2
                node.add_face(support_face, column=0, position="branch-top")
    
    print(f"Applied bold formatting to {formatted_count} placed sequences")
    
    # Set title
    title_face = TextFace(f"{supergroup_of_interest} Reference Tree with pplacer Placed Sequences" if supergroup_of_interest else "Reference Tree with pplacer Placed Sequences", fsize=16, bold=True)
    ts.title.add_face(title_face, column=0)
    
    # Add legend
    if supergroup_of_interest:
        legend_text = f"Tree rooted with {supergroup_of_interest} as ingroup\nBlack/Bold: Placed sequences ({len(placed_sequences)})\nBlack: Reference sequences ({num_taxa - len(placed_sequences)})"
    else:
        legend_text = f"Black/Bold: Placed sequences ({len(placed_sequences)})\nBlack: Reference sequences ({num_taxa - len(placed_sequences)})"
    legend_face = TextFace(legend_text, fsize=10)
    ts.legend.add_face(legend_face, column=0)
    
    # Render and save
    print(f"Rendering tree to: {output_file}")
    tree.render(output_file, w=figure_width, h=figure_height, tree_style=ts, dpi=300)
    
    print(f"Tree successfully rendered with ete3")
    print(f"Figure size: {figure_width}x{figure_height} pixels")
    print(f"Highlighted {len(placed_sequences)} placed sequences in bold")


def plot_tree(tree_file, output_file, db_path=None, supergroup_of_interest=None, reference_package_dir=None, fasta_file=None):
    """
    Plot phylogenetic tree with placed sequences highlighted using ete3.
    
    Args:
        tree_file: Path to input tree file
        output_file: Path to output image file
        db_path: Optional path to SQLite database
        supergroup_of_interest: Optional supergroup name to filter
        reference_package_dir: Optional path to reference package directory
        fasta_file: Optional path to FASTA file containing sequences for length calculation
    """
    if not os.path.exists(tree_file):
        raise FileNotFoundError(f"Tree file not found: {tree_file}")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    try:
        from ete3 import Tree, TreeStyle, NodeStyle, TextFace, AttrFace
        print("Using ete3 for tree visualization...")
        plot_tree_ete3(tree_file, output_file, db_path, supergroup_of_interest, reference_package_dir, fasta_file)
        
    except ImportError:
        raise ImportError("ete3 is required for tree plotting. Please install with: mamba install ete3")
    except Exception as e:
        raise RuntimeError(f"Tree plotting failed: {e}")

def main():
    """Main function for command line usage."""
    if len(sys.argv) < 3:
        print("Usage: python plot_tree.py <tree_file> <output_file> [database_file] [supergroup_of_interest] [reference_package_dir] [fasta_file]")
        print("  tree_file: Path to input tree file (Newick format)")
        print("  output_file: Path to output image file (PDF, PNG, SVG)")
        print("  database_file: Optional path to SQLite database for identifying placed sequences")
        print("  supergroup_of_interest: Optional supergroup name to use for rooting the tree")
        print("  reference_package_dir: Optional path to reference package directory for reference sequence taxonomy")
        print("  fasta_file: Optional path to FASTA file containing sequences for length calculation (for placed sequences only)")
        sys.exit(1)
    
    tree_file = sys.argv[1]
    output_file = sys.argv[2]
    db_path = sys.argv[3] if len(sys.argv) > 3 else None
    supergroup_of_interest = sys.argv[4] if len(sys.argv) > 4 else None
    reference_package_dir = sys.argv[5] if len(sys.argv) > 5 else None
    fasta_file = sys.argv[6] if len(sys.argv) > 6 else None
    
    # If no reference package directory provided, try to find it based on supergroup
    if not reference_package_dir and db_path:
        # Try to find reference package in common locations
        script_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Use supergroup of interest if provided, otherwise default to Amoebozoa
        ref_supergroup = supergroup_of_interest if supergroup_of_interest else 'Amoebozoa'
        possible_refpkg = os.path.join(script_dir, f'../data/reference_packages/{ref_supergroup}.refpkg')
        
        if os.path.exists(possible_refpkg):
            reference_package_dir = possible_refpkg
            print(f"Found reference package at: {reference_package_dir}")
        else:
            # Fallback to Amoebozoa if the specific supergroup package doesn't exist
            fallback_refpkg = os.path.join(script_dir, '../data/reference_packages/Amoebozoa.refpkg')
            if os.path.exists(fallback_refpkg):
                reference_package_dir = fallback_refpkg
                print(f"Supergroup {ref_supergroup} package not found, using fallback: {reference_package_dir}")
    
    if supergroup_of_interest:
        print(f"Rooting tree with {supergroup_of_interest} as ingroup")
    
    plot_tree(tree_file, output_file, db_path, supergroup_of_interest, reference_package_dir, fasta_file)

if __name__ == "__main__":
    main()