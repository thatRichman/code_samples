from chemlinks import get_mechanism
from uuid import uuid4

class DBObj(object):
    """Base class for various database objects (drug, target, associations)

    DBObj and its derivative classes make use of __slots__, greatly decreasing the memory requirements and speed when thousands of these objects
    are instantiated

    meta-tuples of what slots are considered atomic identifiers and what are not are included in the slots 'atomics' and 'non_atomics'
    """
    __slots__ = (
        'uuid',
        '_destroy', # flag used when unifying
        'src', # flag used when unifying
        'parent', # used when unifying
        # atomics
        'atomics',
        'alfred_id',
        'cas_number',
        'chebi_id',
        'chembl_id',
        'comp_tox_id',
        'ensembl_id',
        'gen_atlas_id',
        'hgnc_symbol',
        'humancyc_gene_id',
        'ncbi_gene_id',
        'omim_id',
        'pharmgkb_id',
        'ps_id',
        'pubchem_cid',
        'refseq_dna_id',
        'refseq_prot_id',
        'refseq_rna_id',
        'ttdb_id',
        'uniprot_accession',
        # non-atomics
        'non_atomics',
        'add_ensembl_ids',
        'add_uniprot_accs',
        'aliases',
        'chembl_type',
        'max_phase',
        'organism',
        'pharmgkb_type',
        'pubchem_sids',
        'ttdb_type',
        'has_tract_dat',
        'ttdb_name',
        'pharmgkb_name',
        'chembl_name',
        'ot_name',
        'other_symbols',
    )

    def __init__(self) -> None:
        self.uuid = str(uuid4())
        self._destroy = False
        self.src = None
        self.parent = None
        # store self-references to atomic and non-atomic slot names
        self.atomics = ( # update this when updating __slots__
        'alfred_id',
        'cas_number',
        'chebi_id',
        'chembl_id',
        'comp_tox_id',
        'ensembl_id',
        'gen_atlas_id',
        'hgnc_symbol',
        'humancyc_gene_id',
        'ncbi_gene_id',
        'omim_id',
        'pharmgkb_id',
        'ps_id',
        'pubchem_cid',
        'refseq_dna_id',
        'refseq_prot_id',
        'refseq_rna_id',
        'ttdb_id',
        'uniprot_accession',)
        self.non_atomics = ( # update this when updating __slots__
        'add_ensembl_ids',
        'add_uniprot_accs',
        'aliases',
        'chembl_type',
        'max_phase',
        'organism',
        'pharmgkb_type',
        'pubchem_sids',
        'ttdb_type',
        'has_tract_dat',
        'ttdb_name',
        'pharmgkb_name',
        'chembl_name',
        'ot_name',
        'other_symbols',
        )

    def __iter__(self):
        if not hasattr(self, 'mro'):
            ob = type(self)
        else:
            ob = self
        all_slots = []
        for cls in reversed(ob.mro()): # work backwards in inheritance tree to get all __slots__ and combine
            all_slots.extend(cls.__dict__.get('__slots__', ()))
        tuples = [(a , getattr(self, a, None)) for a in all_slots]
        for t in tuples:
            yield t

    def __eq__(self, o: object) -> bool:
        """Equality operator overidden.

        If the atomic identifier sets of two objects are disjoint, return False. Else return True.

        Args:
            o (object): object to determine equality against. Must define an `atomic_ids_not_null()` method

        Returns:
            bool: whether or not the objects are equal
        """
        if not set(self._atomic_ids_not_null()).isdisjoint(o._atomic_ids_not_null()): # signifcantly faster when tested against generator and intersection
            return True
        return False

    def __repr__(self) -> str:
        return f"Instance of class {type(self).__name__} with {len(list(self._atomic_ids_nn()))} atomic identifiers"

    # tuple of atomic ids
    def _atomic_ids(self):
        return (getattr(self, k, None) for k in self.atomics)

    def _atomic_ids_not_null(self):
        return [x for x in [getattr(self, k, None) for k in self.atomics] if x is not None] # faster than sub-calling atomic_ids

    # tuple of non-atomic values
    def _non_atomics(self):
        return (getattr(self, k, None) for k in self.non_atomics)

    # tuple of non-null non-atomic values
    def _non_atomics_not_null(self):
        return [k for k in self._non_atomics() if k is not None]

    # return all slots of object as dict
    def to_dict(self):
        return {k[0] : k[1] for k in self.__iter__() if k[0] != 'parent'}

    # return all atomic slots of object as dict
    def atomics_dict(self):
        return {k : getattr(self, k, None) for k in self.atomics}

    # return all non-atomic slots of object as dict
    def non_atomics_dict(self):
        return {k : getattr(self, k, None) for k in self.non_atomics}

    def node_dict(self):
        return {k[0] : k[1] for k in self.__iter__() if k[0] not in ['parent','_destroy','uuid','atomics','non_atomics','_']}

    def to_node(self):
        return (self.uuid, self.node_dict())

    def merge(self, y):
        """merge object y into self

        Currently included for posterity
        """
        if type(self).__name__ != type(y).__name__:
            raise ValueError("These two objects are not the same class or subclass")
        for k, v in y.atomics_dict().items():
            if v:
                setattr(self, k, v)
        for k,yattr in y.non_atomics_dict().items():
            if yattr is not None:
                xattr = getattr(self, k, None)
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
                    setattr(self, k, yattr)

        if isinstance(self.src, list):
            if isinstance(y.src, list):
                self.src.extend(y.src)
            else:
                self.src.append(y.src)
        else:
            if isinstance(y._src, list):
                y._src.append(self.src)
                self.src = y._src
            else:
                self.src = [self.src, y.src]
        self.uuid = str(uuid4())


class Drug(DBObj):
    """Drug subclass of DBObj
    """
    __slots__ = ('_',)

    def __init__(self) -> None:
        super().__init__()

    @classmethod
    def from_dict(self, dct):
        d = Drug()
        for k, v in dct.items():
            try:
                setattr(d, k, v)
            except AttributeError as ae:
                pass
        return d


class Target(DBObj):
    """Target subclass of DBObj
    """
    __slots__ = ('_',)

    def __init__(self) -> None:
        super().__init__()

    @classmethod
    def from_dict(self, dct):
        t = Target()
        for k, v in dct.items():
            setattr(t, k, v)
        return t

    def __uniprot_eq__(self, o: object) -> bool:
        if set(getattr(self, 'add_uniprot_ids',[])) & set(getattr(o, 'add_uniprot_ids', [])):
            return True
        return False

    def __ensembl_eq__(self, o: object) -> bool:
        if set(getattr(self, 'add_ensembl_accs',[])) & set(getattr(o, 'add_ensembl_accs',[])):
            return True
        return False


class Association(DBObj):
    """Association subclass of DBObj
    """
    __slots__ = (
        'drug_uuid',
        'target_uuid',
        'mechanism',
        'mechanism_text',
        'mechanism_original',
        'direct_interaction',
        'conflict_resolved',
    )

    def __repr__(self) -> str:
        return f"Association object between drug uuid {self.drug_uuid} and target uuid {self.target_uuid}"

    def __init__(self, drug_uuid, target_uuid, **kwargs) -> None:
        self.uuid = str(uuid4())
        self.drug_uuid = drug_uuid
        self.target_uuid = target_uuid
        self.mechanism = kwargs.get('mechanism')
        self.conflict_resolved = None
        for k, v in kwargs.items():
            try:
                setattr(self, k, v)
            except AttributeError:
                pass

    def __eq__(self, o: object) -> bool:
        """Special case of DBOBj equality

        Returns True if and only if both the drug uuid of self and object, and target uuid of self and object are equal

        Args:
            o (object): object to be tested for equality against

        Returns:
            bool: whether the two objects are equal
        """
        if self.drug_uuid == o.drug_uuid and self.target_uuid == o.target_uuid:
            return True
        return False


