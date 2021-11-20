import json
import multiprocessing
from functools import partial
from itertools import combinations, islice
from multiprocessing.queues import Queue
from time import process_time
from uuid import uuid4

import networkx as nx
from networkx.classes.function import edges
from networkx.readwrite import json_graph
from numpy.core.function_base import add_newdoc
from tqdm import tqdm

from conflict_resolution import (consensus_resolve, hierarchy_resolve,
                                 proportional_resolve)


class IdGraph:
    """Graph representation for reconciliation of drugs or targets
    """
    def __init__(self, dbobjs = None) -> None:
        self.graph = nx.Graph()
        self.uuids = []
        self.node_dict = {}
        dbobjs = [dbo for dbo in dbobjs] # consume the generator
        self.dbobjs = dbobjs
        if dbobjs:
            for obj in dbobjs:
                self.node_dict[obj.uuid] = obj
            self.add_nodes(dbobjs)

    def add_nodes(self, dbobjs) -> None:
        """Add DBObjects to the graph

        Args:
            dbobjs (DBObj or list): iterable of dbobjs
        """
        nodes = (obj.to_node() for obj in dbobjs)
        for uuid, node in tqdm(nodes):
            if uuid not in self.uuids:
                self.graph.add_node(node_for_adding=uuid, **node)
                self.uuids.append(uuid)

    def auto_edges(self) -> None:
        """
        Automatically add edges between existing graph nodes based on identifier set disjoint
        """
        # fastest in practice. Trotter's recursive algorithm is too much overhead, and pre-listing is often too large
        # combinations are the speed bottleneck
        # using integer indices uses less memory but isn't any faster
        uuids = tuple(self.graph.nodes)
        combos = combinations(range(len(uuids)),2)
        edges_to_add = []
        for x in tqdm(combos):
            if self.node_dict[uuids[x[0]]] == self.node_dict[uuids[x[1]]]:
                edges_to_add.append((x[0],x[1]))
        self.graph.add_edges_from(edges_to_add)

    def make_edges(self, nodes) -> None:
        """Check if any edges should exist between nodes, and if so, make them

        Args:
            nodes ([type]): [description]
        """
        for x in tqdm(combinations(nodes, 2)):
            if self.node_dict[x[0]] == self.node_dict[x[1]]:
                self.graph.add_edge(x[0], x[1])

    def edges_from_tuples(self, tuples, indices = False):
        """Check if an edge should exist between elements in each 2-tuple, and if so, make them

        Args:
            tuples (iterable): iterable of tuples
            indices (bool): if True, tuple values are treated as indices of self.graph.nodes, else, treated as node names (default)
        """
        print(f'checking {len(tuples)} potential edges')
        if not indices:
            for x in tqdm(tuples):
                if x[0] == x[1]:
                    self.graph.add_edge(x[0],x[1])
        else:
            nodes = tuple(self.graph.nodes)
            for x in tqdm(tuples):
                if  self.node_dict[nodes[x[0]]] == self.node_dict[nodes[x[1]]]:
                    self.graph.add_edge(nodes[x[0]], nodes[x[1]])

    def merge_connected(self):
        """Contract all maximal groups of connected nodes into single nodes

        Returns:
            list: contracted nodes
            dict: mapping between input node uuids and contracted node uuids
        """
        connected = nx.connected_components(self.graph)
        for c in tqdm(connected):
            to_merge = (self.node_dict[x] for x in c)
            md = to_merge[0]
            if len(to_merge) > 1:
                for d in islice(to_merge, 1):
                    md = self.merge(md, d)
            yield md

    def merge(self, x, y):
        """combine two objects into one"""
        if type(x).__name__ != type(y).__name__:
            raise ValueError("These two objects are not the same class or subclass")
        for k, v in y.atomics_dict().items():
            if v:
                setattr(x, k, v)
        for k,yattr in y.non_atomics_dict().items():
            if yattr is not None:
                xattr = getattr(x, k, None)
                if xattr is not None:
                    if isinstance(xattr, list) and isinstance(yattr, list):
                        xattr.extend(yattr)
                    elif isinstance(xattr, list):
                        xattr.append(yattr)
                    elif isinstance(yattr, list):
                        yattr.append(xattr)
                        xattr = yattr.copy()
                    else:
                        xattr = [xattr, yattr]
                else:
                    setattr(x, k, yattr)

        if isinstance(x.src, list):
            if isinstance(y.src, list):
                x.src.extend(y.src)
            else:
                x.src.append(y.src)
        else:
            if isinstance(y.src, list):
                y.src.append(x)
                x.src = y.src
            else:
                x.src = [x.src, y.src]
        x.uuid = str(uuid4()) # regenerate new uuid
        return x

    def update(self, objs):
        """Update graph with new edges
        Args:
            objs (DBObj): New Objects to be added
        """
        self.add_nodes(objs)
        isos = nx.isolates(self.graph)
        # check all isolates against each other
        self.make_edges(isos)
        # check isolates against connected
        connected_nodes = (x[0] for x in self.graph.degree if x[1] >= 1 and x[0] not in isos)
        edges_to_check = [(x,y) for x in isos for y in connected_nodes]
        self.edges_from_tuples(edges_to_check)


class AssocGraph(IdGraph):
    """Subclass of IdGraph specifically for working with Association objects
    """
    def __init__(self, dbobjs) -> None:
        super().__init__(dbobjs=dbobjs)

    # associations require a conflict resolution layer
    def merge_connected(self, check_conflicts=True, resolve_method='hybrid'):
        """Contract maximal groups of connected Associations based on drug / target uuids

        Returns:
            list: contracted Associations
            dict: mapping of input Associations to contracted Associations
        """
        connected = nx.connected_components(self.graph)
        for c in tqdm(connected):
            if len(c) > 1:
                # is there a conflict?
                to_merge = [self.node_dict[x] for x in c]
                if check_conflicts:
                    mechs = [getattr(x, 'mechanism') for x in to_merge]
                    if mechs.count(mechs[0]) != len(mechs):  # conflict
                        to_merge = self.resolve_conflicts(to_merge, resolve_method=resolve_method)
                        mn.conflict_resolved = True
                mn = to_merge[0]
                for n in to_merge[1:]:
                    mn = self.merge(mn, n)
                yield(mn)
            else:
                yield self.node_dict[list(c)[0]]

    # merging associations is a slightly different process than drugs / targets
    def merge(self, x, y):
        """combine two associations into one"""
        if type(x).__name__ != type(y).__name__:
            raise ValueError("These two objects are not the same class or subclass")
        for k, yattr in y.to_dict().items():
            if yattr is not None:
                xattr = getattr(x, k, None)
                if xattr is not None:
                    if isinstance(xattr, list) and isinstance(yattr, list):
                        xattr.extend(yattr)
                    elif isinstance(xattr, list):
                        xattr.append(yattr)
                    elif isinstance(yattr, list):
                        yattr.append(xattr)
                        xattr = yattr.copy()
                    else:
                        xattr = [xattr, yattr]
                else:
                    setattr(x, k, yattr)

        if isinstance(x.src, list):
            if isinstance(y.src, list):
                x.src.extend(y.src)
            else:
                x.src.append(y.src)
        else:
            if isinstance(y._src, list):
                y.src.append(x)
                x.src = y.src
            else:
                x.src = [x.src, y.src]
        x.uuid = str(uuid4()) # regenerate new uuid
        return x

    def resolve_conflicts(self, nodes, resolve_method='hybrid'):
        """Resolve conflicts between associations using specified method

        Currently supports 'hierarchy','consensus', and 'hybrid' / 'proportional'

        Associations with incorrect mechanisms (as determined by method) are culled

        Args:
            nodes (list): list of Associations with conflicts
            method (str, optional): conflict resolution method. Defaults to 'hybrid'.

        Raises:
            ValueError: raised if invalid method provided

        Returns:
            list: list of nodes with correct mechanism as determined by selected method
        """
        if resolve_method == 'consensus':
            return consensus_resolve(nodes)
        elif resolve_method == 'hierarchy':  # arbitrary hierarchy defines which source is more correct
            return hierarchy_resolve(nodes)
        elif resolve_method == 'hybrid' or resolve_method == 'proportional':
            return proportional_resolve(nodes)
        else:
            raise ValueError('Unrecognized conflict resolution method')


class DTGraph:
    """Graph representation of drug-target associations
    """
    def __init__(self, drugs=None, targets=None, associations=None) -> None:
        self.graph = nx.DiGraph()
        self.drugs = {}
        self.targets = {}
        self.associations = {}
        if drugs:
            for drug in drugs:
                self.drugs[drug.uuid] = drug
                self.graph.add_node(node_for_adding=drug.uuid, **drug.to_dict())
        if targets:
            for target in targets:
                self.targets[target.uuid] = target
                self.graph.add_node(node_for_adding=target.uuid, **target.to_dict())
        if associations:
            for association in associations:  # associations are not added as nodes
                self.associations[association.uuid] = association

    def make_edges(self):
        """Link drugs to targets based on associations
        """
        if not self.drugs or not self.targets:
            print("no edges to make")
            return
        elif not self.associations:
            print("no associations from which to make edges")
            return
        for _, a in self.associations.items():
            if a.drug_uuid in self.graph.nodes and a.target_uuid in self.graph.nodes:
                self.graph.add_edge(a.drug_uuid, a.target_uuid, **a.to_dict())  # drug --> target

    def to_json(self, path):
        """Export graph as JSON file

        Args:
            path (str or Path): Location for writing
        """
        jdat = json_graph.node_link_data(self.graph)
        with open(path, 'w') as jf:
            jf.writelines(json.dumps(jdat))

    def to_json_cytoscape(self, path=None):
        """Export graph as Cytoscape-compatible JSON file

        Args:
            path (str or Path): Location for writing
        """
        jdat = json_graph.cytoscape.cytoscape_data(self.graph)
        if path:
            with open(path, 'w') as jf:
                jf.writelines(json.dumps(jdat))
        return jdat
