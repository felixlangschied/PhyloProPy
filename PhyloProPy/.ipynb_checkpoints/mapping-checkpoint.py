
def check_taxonomy_input(lineage, ncbi):
    """Returns taxid as int"""
    # check format of lineage
    if isinstance(lineage, int):
        return lineage
    elif str(lineage).startswith('ncbi'):
        return int(lineage.replace('ncbi', ''))
    else:
        # try to parse species name to taxid
        name2taxid = ncbi.get_name_translator([lineage, lineage.replace('_', ' ')])
        if len(name2taxid) == 1:
            return list(name2taxid.values())[0][0]