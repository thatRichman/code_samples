import argparse
import importlib
import time
import warnings
from pathlib import Path
from typing import Type
import pandas as pd

from gooey import Gooey, GooeyParser

from chembl import chembl
from DBObj import Association, Drug, Target
from id_graph import AssocGraph, IdGraph, DTGraph
from opentargets import opentargets
from pharmgkb import pharmgkb
from ttdb import ttdb
from DrugDB import DrugDB

from tqdm import tqdm

SOURCES = ('chembl', 'ttdb', 'opentargets', 'pharmgkb') # default sources, don't necessarily need to append new sources here


def preprocess(sources, check_files=True):
    """Execute the `run` method for each Source

    `run` should be declared for every Source, and execute the necessary processing steps to turn
    raw data from the Source into correctly-formatted JSON files.

    Args:
        sources (list): list of sources (as strs) to be used in database construction

    Returns:
        list: Source object for each name in Source
    """
    objs = [eval(source)() for source in sources]

    for o in tqdm(objs):
        o.run(check_files = check_files)
    return objs


# def make_dbobjs(objs):
#     """Convert dicts loaded from JSON into Database Objects

#     Args:
#         objs (list): list of dicts loaded from JSON file

#     Returns:
#         dict: Dict containing 'drugs', 'targets', 'associations' entries, each of which is a list of DBObjs of appropriate subclass
#     """
#     drugs = []
#     targets = []
#     associations = []
#     for o in objs:
#         src = type(o).__name__
#         print(f'Loading {src} objects')
#         d_tmp = o.load_drugs()
#         # Technically, all sources will have all of these methods defined. So cannot use hasattr as a test. Instead, just check if return value is not None
#         if d_tmp:
#             d_objs = [Drug.from_dict(d) for d in d_tmp]
#             for d in d_objs:
#                 d.src = src
#             drugs.extend(d_objs)
#         t_tmp = o.load_targets()
#         if t_tmp:
#             t_objs = [Target.from_dict(t) for t in t_tmp]
#             for t in t_objs:
#                 t.src = src
#             targets.extend(t_objs)
#         a_tmp = o.load_associations()  # Association dbobjs are not built here, dicts are just loaded
#         if a_tmp:
#             for a in a_tmp:
#                 a['src'] = src
#             associations.extend(a_tmp)
#     return {
#         'drugs' : drugs,
#         'targets' : targets,
#         'associations' : associations
#     }

def make_drugs(objs):
    for o in tqdm(objs):
        src = type(o).__name__
        d_tmp = o.load_drugs()
        if d_tmp:
            d_objs = [Drug.from_dict(d) for d in d_tmp]
            for d in d_objs:
                d.src = src
                yield d


def make_targets(objs):
    for o in tqdm(objs):
        src = type(o).__name__
        t_tmp = o.load_targets()
        if t_tmp:
            t_objs = [Target.from_dict(t) for t in t_tmp]
            for t in t_objs:
                t.src = src
                yield t


def make_associations(objs):
    for o in tqdm(objs):
        src = type(o).__name__
        a_tmp = o.load_associations()  # Association dbobjs are not built here, dicts are just loaded
        if a_tmp:
            for a in a_tmp:
                a['src'] = src
                yield a


def unify(db_objs):
    """Reconcile drugs or targets from varying sources

    Args:
        db_objs (iterable): iterable of database objects (of same subclass) to be unified

    Returns:
        unified (list): list of unified drugs
        uuid_dict (dict): dict mapping original uuids of input to uuids of output
    """
    igraph = IdGraph(db_objs)
    igraph.auto_edges()
    return igraph.merge_connected()


def build_associations(associations, drugs, targets):
    """Convert association dictionaries in Association objects with correct drug / target uuid mapping

    Associations that do not have both a matching drug and target uuid are skipped
    Associations whose src is unsupported are also skipped

    Args:
        associations (list): list of dicts to be converted
        drugs (list): drug objects from which to search for required uuids
        targets (list): target objects from which to search for required uuids

    Returns:
        assoc_objs (list): list of Association Objects
    """
    for assoc in tqdm(associations):
        if assoc['src'] == 'chembl':
            drug_uuid = [d.uuid for d in drugs if getattr(d, 'chembl_id', None) == assoc['drug_chembl_id']]
            target_uuid = [t.uuid for t in targets if getattr(t, 'chembl_id', None) == assoc['target_chembl_id']]
            if drug_uuid and target_uuid:
                drug_uuid = drug_uuid[0]
                target_uuid = target_uuid[0]
                yield Association(drug_uuid = drug_uuid, target_uuid = target_uuid, **assoc)
        else:
            raise ValueError(f"Unsupported source {assoc['src']}, skipping")


def unify_associations(associations, check_conflicts = True, resolve_method='hybrid'):
    """Reconcile associations from various sources and resolve conflicts

    Args:
        associations (list): list of Association objects to be unfied

    Returns:
        unified (list): list of unified Association objects
        uuid_dict (dict): dict mapping original uuids to unified uuids
    """
    print(f"Unifying {len(associations)} associations")
    agraph = AssocGraph(associations)
    agraph.auto_edges()
    return agraph.merge_connected(check_conflicts = check_conflicts, resolve_method=resolve_method)


def build_dtgraph(drugs, targets, associations):
    """Create Drug-Target association digraph from drug, target, and association objects

    Args:
        drugs (Drug): Drug objects
        targets (Target): Target objects
        associations (Association): Association objects
    """
    dtgraph = DTGraph(drugs = drugs, targets = targets, associations = associations)
    dtgraph.make_edges()
    return dtgraph


def make_sqlite_db(drugs, targets, associations, path):
    """Make sqlite database from drugs, targets, and associations

    Args:
        drugs (Drug): Drug objects
        targets (Target): Target objects
        associations (Association): Association objects
    """

    # Convert objects to dicts
    drug_dicts = [d.to_dict() for d in drugs]
    target_dicts = [t.to_dict() for t in targets]
    assoc_dicts = [a.to_dict() for a in associations]

    # Convert everything to strings
    for idx, val in enumerate(drug_dicts):
        for k,v in val.items():
            if isinstance(v, (list,tuple)):
                drug_dicts[idx][k] = ','.join(v)
    for idx, val in enumerate(target_dicts):
        for k,v in val.items():
            if isinstance(v, (list,tuple)):
                target_dicts[idx][k] = ','.join(v)
    for idx, val in enumerate(assoc_dicts):
        for k,v in val.items():
            if isinstance(v, (list,tuple)):
                assoc_dicts[idx][k] = ','.join(v)

    db = DrugDB(path)
    db.build(drug_dicts, target_dicts, assoc_dicts)
    return


def nrdt_args(subparsers):
    """Run the GUI."""
    argparser = subparsers.add_parser(
        'NRDT_args',
        prog='NRDT Database Construction',
        description="""
        Build Non-redundant drug-target database.
        This is a very slow process and should only be run if not provided a preconstructed database.
        """,
        )
    argparser.add_argument(
        'sources',
        help='Select sources to include in database build',
        widget='Listbox',
        choices=SOURCES,
        nargs='+',
        gooey_options = {
            'height': 100}

    )
    argparser.add_argument(
        'path',
        help = 'Path to write SQLite DB',
    )
    argparser.add_argument(
        '--resolve_method',
        help = 'Choose conflict resolution method',
        choices = ['proportional','consensus','hierarchy'],
        default = 'proportional'
    )
    argparser.add_argument(
        '--check_files',
        help = 'Check each source for existing fully-processed files. If exist, will not reprocess',
        action='store_true',
        default=True
    )
    argparser.add_argument(
        '--check_conflicts',
        help = 'Check for conflicts in interactions',
        action = 'store_true',
        default = True
    )
    # argparser.add_argument(
    #     '--test',
    #     help = 'Internal testing. If you see this option, bug the developer',
    #     action = 'store_true',
    #     default = False
    # )
    argparser.set_defaults(func=build_db)


def build_db(args):
    pd.options.mode.chained_assignment = None  # default='warn'
    db_path = Path(args.path)
    for source in args.sources:
        try:
            importlib.import_module(source, package=source)
        except ImportError:
            print(f'Failed to load source {source}')
            exit(1)
    print('Preprocessing')
    start = time.perf_counter()
    src_procd = preprocess(args.sources, args.check_files)
    stop = time.perf_counter()
    print(f'Preprocessed {len(SOURCES)} sources in {stop-start} seconds')
    print('Making objects')
    start = time.perf_counter()
    drugs = make_drugs(src_procd)
    targets = make_targets(src_procd)
    associations = make_associations(src_procd)
    stop = time.perf_counter()
    # reconcile drugs and targets
    print("Unifying drugs...")
    uni_drugs = unify(drugs)
    print("Unifying targets...")
    uni_targs = unify(targets)
    # build and reconcile associations
    print("Unifying associations...")
    associations = [a for a in build_associations(associations, drugs, targets)]
    if not args.check_conflicts:
        warnings.warn("Not checking for conflicts")
    uni_assocs = unify_associations(associations, check_conflicts = args.check_conflicts, resolve_method = args.resolve_method)
    # Build Drug-Target DiGraph
    print("Building DTGraph...")
    dtgraph = build_dtgraph(uni_drugs, uni_targs, uni_assocs)
    # dump DTGraph to json
    print("Dumping DTGraph to JSON...")
    jdat = dtgraph.to_json('dtgraph.json')
    # convert to sqlite
    print("Building SQlite database...")
    make_sqlite_db(drugs = uni_drugs, targets = uni_targs, associations = uni_assocs, path = db_path)
    print("...done.")

@Gooey(
    program_name='NRDT-DB Builder',
    advanced=True,
    default_size=(1024, 768),
    navigation='TABBED',
    show_stop_warning=False,
    show_restart_button=True,
    clear_before_run=True,
    progress_regex=r'(?P<current>\d+) / (?P<total>\d+)$',
    progress_expr='current / total * 100',
    hide_progress_msg=True,
)
def main():
    argparser = GooeyParser(
        description='Construct Non-Redundant Drug-Target Interaction Database',
    )
    subparsers = argparser.add_subparsers()
    nrdt_args(subparsers)

    for choice in subparsers.choices.values():
        for action in choice._actions:
            if action.help:
                action.help = ' '.join(action.help.split())

    args = argparser.parse_args()
    if not args.sources:
        args.sources = ["ttdb","opentargets","chembl","pharmgkb"]
    args.func(args)


if __name__ == '__main__':
    main()
