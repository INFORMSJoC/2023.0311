import os
import json
import sys

import networkx as nx
import numpy as np
import pandas as pd

from heapq import heappush, heappop
from itertools import count

PATH_COLS = ["path_id", "departure_min", "cost"]
GRAPH_COLS = ["origin_trip_id", "dest_trip_id", "edge_cost"]

SOURCE_ID = "src"
SINK_ID = "sink"
DROPOFF_SUFFIX = "_dropoff"
PICKUP_SUFFIX = "_pickup"



config_filepath = sys.argv[1]
with open(config_filepath, 'r') as file: 
    config = json.load(file)

OUTPUT_DIR = config["output_dir"]
MAX_PATHS = config["num_itin"]
HOUR_WAGE = config["hourly_wage"]
SEC_PER_MIN = 60
MINUTE_WAGE = HOUR_WAGE / SEC_PER_MIN

####################### GENERATE ITINERARIES

"""
Yen's k shortest loopless paths algorithm. Requires non-negative edge weights.

Adapted from: https://github.com/guilhermemm/k-shortest-path/blob/master/k_shortest_paths.py

### Keywords
* `G`
* `source`
* `target`
* `trip_nodes`
* `k`
* `weight`

### Returns
* path lengths
* paths
"""
def k_shortest_paths(G, source, target, k=1, weight='cost'):
    if source == target:
        return ([0], [[source]])

    try:
        length, path = nx.single_source_dijkstra(G, source, target, weight=weight)
    except nx.NetworkXNoPath:
        print("node %s not reachable from %s" % (source, target))
        raise

    paths = [path]
    c = count()
    B = []
    G_original = G.copy()

    for i in range(1, k): # 1 to k-1
        # spur node ranges from the first node to the next to last node in the previous k-shortest path.
        for j in range(len(paths[-1]) - 1): # for each root path
            spur_node = paths[-1][j]
            root_path = paths[-1][:j + 1]

            edges_removed = []
            for c_path in paths:
                # remove the next-links that are part of the previous shortest paths which share the same root path
                if len(c_path) > j and root_path == c_path[:j + 1]:
                    u = c_path[j]
                    v = c_path[j + 1]
                    if G.has_edge(u, v):
                        edge_attr = G.edges[u,v]
                        G.remove_edge(u, v)
                        edges_removed.append((u, v, edge_attr))

            # disconnect every root path node except spur node from the rest of the graph to ensure simple path
            for n in range(len(root_path) - 1):
                node = root_path[n]

                out_edges=list(G.out_edges(node))
                for u, v in out_edges:
                    edge_attr=G.edges[(u,v)]
                    G.remove_edge(u, v)
                    edges_removed.append((u, v, edge_attr))

                if G.is_directed():
                    in_edges = list(G.in_edges(node))
                    for u, v in in_edges:
                        edge_attr=G.edges[(u,v)]
                        G.remove_edge(u, v)
                        edges_removed.append((u, v, edge_attr))

            # find shortest spur path, guaranteed to be simple
            try:
                spur_path_length, spur_path = nx.single_source_dijkstra(G, spur_node, target, weight=weight)
            except nx.NetworkXNoPath:
                # restore graph, move on
                for e in edges_removed:
                    u, v, edge_attr = e
                    G.add_edge(u, v)
                    G[u][v].update(edge_attr)
                continue
            else:
                # record full path
                total_path = root_path[:-1] + spur_path
                total_path_length = get_path_length(G_original, root_path, weight) + spur_path_length
                heappush(B, (total_path_length, next(c), total_path))


            # restore graph
            for e in edges_removed:
                u, v, edge_attr = e
                G.add_edge(u, v)
                G[u][v].update(edge_attr)

        if B:
            (_, _, p) = heappop(B)
            # take care of duplicate paths
            while p in paths and B:
                (_, _, p) = heappop(B)
            if p not in paths:
                paths.append(p)
        else:
            break

    return paths


def get_path_length(G, path, weight='cost'):
    length = 0
    if len(path) > 1:
        for i in range(len(path) - 1):
            u = path[i]
            v = path[i + 1]

            length += G[u][v].get(weight, 1)

    return length

"""
Read in network.

### Keyword Args
* `output_dir` - directory containing all files
* `graph_fn` - filename of graph file
* `trip_fn` - filename of trip file
### Returns
* NetworkX DiGraph
* trip DataFrame
"""
def itinerary_init(output_dir, graph_fn, trip_fn):
    trip = pd.read_csv(os.path.join(output_dir, trip_fn))
    graph_dat = pd.read_csv(os.path.join(output_dir, graph_fn))
    assert all(c in graph_dat.columns for c in GRAPH_COLS)

    # build graph
    G = nx.DiGraph()

    #--- expand trip nodes
    for i in range(len(trip)):
        row = trip.iloc[i]
        origin = str(row['trip_id']) + PICKUP_SUFFIX
        dest = str(row['trip_id']) + DROPOFF_SUFFIX
        G.add_edge(origin, dest, cost=row['travel_time_min'] * MINUTE_WAGE)

    #--- only add edges meant to build itineraries for planned drivers
    for i in range(len(graph_dat)):
        row = graph_dat.iloc[i]
        if row['planned_driver_edge']:
            origin = SOURCE_ID if row['origin_trip_id'] == SOURCE_ID else str(row['origin_trip_id']) + DROPOFF_SUFFIX
            dest = SINK_ID if row['dest_trip_id'] == SINK_ID else str(row['dest_trip_id']) + PICKUP_SUFFIX
            G.add_edge(origin, dest, cost=row['edge_cost'])

    return G, trip

"""
Given k shortest paths:
* remove duplicates
* record
"""
def postprocess_ksp(trip, paths, fp):
    # trip IDs
    trip_id = [str(t) for t in trip.trip_id]
    n = len(trip_id)

    # seq of trip pickups
    path_ids = []
    for p in paths:
        p = [n.replace(PICKUP_SUFFIX, '') for n in p]
        p = [n.replace(DROPOFF_SUFFIX, '') for n in p]
        new_p = []
        for n in p:
            if n not in new_p:
                new_p.append(n)
        path_ids.append('-'.join([n for n in new_p]))

    # itinerary candidates - get unique
    path_df = pd.DataFrame(path_ids, columns=['path_id'])

    #--- ensure all trips represented
    unrepresented = [t for t in trip_id if not any([t in path_id for path_id in path_ids])]
    new_paths = []
    new_paths = np.array(['-'.join([SOURCE_ID, n, SINK_ID]) for n in unrepresented]).reshape(-1,1)
    new_path_df = pd.DataFrame(new_paths, columns=['path_id'])
    path_df = pd.concat([path_df, new_path_df]).reset_index(drop=True)
    path_df.to_csv(fp, index=False)
    return


def main():
    trip_fn = 'trips.csv'
    graph_fn = 'edges.csv'
    output_fn = 'itin.csv'
    output_fp = os.path.join(OUTPUT_DIR, output_fn)
    G, trip = itinerary_init(OUTPUT_DIR, graph_fn, trip_fn)

    print("k shortest paths...")
    paths = k_shortest_paths(G, SOURCE_ID, SINK_ID, MAX_PATHS)

    print("post-processing paths...")
    postprocess_ksp(trip, paths, output_fp)

    print("done with k paths.")


if __name__ == '__main__':
    main()
