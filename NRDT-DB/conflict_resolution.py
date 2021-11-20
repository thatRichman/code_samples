from collections import Counter
import warnings

HIERARCHY_DICT = {
    'chembl' : 0,
    'opentargets' : 1,
    'pharmgkb' : 2,
    'ttdb' : 3,
}
PROPORTION_DICT = {
    'chembl' : 4,
    'opentargets' : 3,
    'pharmgkb' : 2,
    'ttdb' : 1,
}


def proportional_resolve(nodes):
    """Perform proportional representation conflict resolution based on source weights in PROPORTION_DICT

    Args:
        nodes (Association): Association nodes with conflicts to be resolved

    Returns:
        list: All nodes with correct mechanism as decided by proportional vote
    """
    out = []
    scores = {}
    for n in nodes:
        weight = PROPORTION_DICT.get(n)
        mech = n.mechanism
        if scores[mech]:
            scores[mech] += weight
        else:
            scores[mech] = weight
    lead = max(scores, key=scores.get)
    for n in nodes:
        if n.mechanism == lead:
            out.append(n)
    return out


def consensus_resolve(nodes):
    """Perform consensus / majority rules conflict resolution

    Args:
        nodes (Association): Association nodes with conflicts to be resolved

    Returns:
        list: All nodes with correct mechanism as decided by consensus / majority rules
    """
    ctr = Counter([n.mechanism for n in nodes])
    top = ctr.most_common(2)
    if top[0][1] == top[1][1]:
        # tie, revert to hierarchy
        warnings.warn('Tie detected, reverting to hierarchical conflict resolution for this association')
        return hierarchy_resolve(nodes)
    else:
        out = []
        for node in nodes:
            if node.mechanism == top[0][0]:
                out.append(node)
        return out


def hierarchy_resolve(nodes):
    """Perform hierarchical conflict resolution based on source ranking in HIERARCHY_DICT

    Args:
        nodes (Association): Association nodes with conflicts to be resolved

    Returns:
        list: All nodes with correct mechanism as decided by hierarchy
    """
    lead = None
    out = []
    h = 100
    for n in nodes:
        nh = HIERARCHY_DICT.get(n.src, 100)
        if nh < h:
            h = nh
            lead = n.mechanism
    for node in nodes:
        if node.mechanism == lead:
            out.append(node)
    return out