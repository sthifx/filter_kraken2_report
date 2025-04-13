#!/usr/bin/env python

"""
filter_kraken2_report
Author: Sam Haldenby
Created: 2025-04-13

Description:
    A Python tool for filtering Kraken2 reports by removing specified taxa and their descendants.
    Uses NetworkX to build and prune the taxonomic tree, then re-exports a trimmed report.

License: GNU General Public License v3.0 (GPL-3.0)
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see https://www.gnu.org/licenses/.

Repository: https://github.com/yourusername/filter_kraken2_report
"""


import argparse
import pandas as pd
import networkx as nx

def remove_tax_and_below(tax_id, g):
    if tax_id not in g.nodes:
        print (f'ID {tax_id} not found in taxonomy tree')
        return g

    total_counts_before = sum([g.nodes[node]['at_rank'] for node in g.nodes])
    target_name = g.nodes[tax_id]['name']

    # Get upstream and downstream nodes
    curr_node = g.nodes[tax_id]
    downstream_nodes = nx.descendants(g, tax_id)
    upstream_nodes = nx.ancestors(g, tax_id)

    # Get at_rank sum values for this node and descendants
    at_node_or_lower = sum([g.nodes[node]['at_rank'] for node in downstream_nodes]) + curr_node['at_rank']

    # Subtract that from all ancestor nodes at_rank_or_lower values
    for node in upstream_nodes:
        g.nodes[node]['at_rank_or_lower'] = g.nodes[node]['at_rank_or_lower'] - at_node_or_lower

    # Remove node and descendents
    downstream_nodes_and_curr = downstream_nodes.copy()
    downstream_nodes_and_curr.add(tax_id)
    g.remove_nodes_from(downstream_nodes_and_curr)

    # Recalculate percentages
    total_counts = sum([g.nodes[node]['at_rank'] for node in g.nodes])
    for node in g.nodes:
        at_rank_or_lower = g.nodes[node]['at_rank_or_lower']
        new_percent = at_rank_or_lower / total_counts * 100.0
        g.nodes[node]['percent'] = new_percent

    # Report
    removed_count = total_counts_before - total_counts
    removed_percentage = removed_count / total_counts_before * 100.0
    print (f'Pruned ID {tax_id} ({target_name}) and children, comprising {removed_count} reads ({removed_percentage :.2f}%)')

    # Return                                                                                                                                                
    return g

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--report", help="kraken2 report file", type=str)
parser.add_argument("-f", "--filter_ids", help="list of tax ids to filter. Comma-separated. These taxa and all directly below will be removed", type=str)
parser.add_argument("-o", "--out_fn", help="output report file name", type=str)
parser.add_argument("-p", "--removed_out_fn", help="output report file name, for recording removed rows", type=str)
opts = parser.parse_args()

# Read report
df = pd.read_csv(opts.report,
                 skipinitialspace=False,
                 sep='\t', 
                 names=['percent','at_rank_or_lower','at_rank','tax_rank','tax_id','name'])

# Calculate the number of leading spaces for each name
df['spaces'] = df['name'].str.extract(r'^( *)')[0].str.len()
df['srank'] = df['spaces'].div(2).astype(int)

# Get tax_id order: This will be used later to reconstruct the output file order
tax_id_order = df['tax_id']

# Create list to contain information on the last parent at each srank/spaces
last_srank = [None] * (df.srank.max()+1)

# Create a graph containing taxa
g = nx.DiGraph()

for _,row in df.iterrows():
    # Set last_srank
    last_srank[row.srank] = row.tax_id
    
    # Add node to graph
    g.add_node(row.tax_id,
              percent = row.percent,
              at_rank_or_lower = row.at_rank_or_lower,
              at_rank = row.at_rank,
              tax_rank = row.tax_rank,
              name = row['name'].strip(),
              raw_name = row['name'])

    # Connect to parent
    if row.srank >0:
        g.add_edge(last_srank[row.srank - 1],row.tax_id)

# Prune
tax_ids = [int(x) for x in opts.filter_ids.split(',')]
for tax_id in tax_ids:
    g = remove_tax_and_below(tax_id, g)
    

# Reconstruct df from graph
rows = []
pruned_ids = []
for tax_id in tax_id_order:
    if tax_id in g.nodes:
        node = g.nodes[tax_id]
        d = dict(percent = f'{node['percent'] :.2f}',
              at_rank_or_lower = node['at_rank_or_lower'],
              at_rank = node['at_rank'],
              tax_rank = node['tax_rank'],
              tax_id = tax_id,
              name = node['raw_name'])
        rows.append(d)
    else:
        pruned_ids.append(tax_id)
    
# Create df
pruned_df = pd.DataFrame(rows)

# Export
pruned_df.to_csv(opts.out_fn, sep='\t', header=False, index=False)

# Export removed entries
removed_df = df.loc[df['tax_id'].isin(pruned_ids),].drop(['spaces','srank'], axis=1)
removed_df.to_csv(opts.removed_out_fn, sep='\t', header=False, index=False)





