from dpgen.generator.arginfo import (
    run_mdata_arginfo,
)
from dpgen.simplify.arginfo import (
    simplify_mdata_arginfo,
)
from dpgen.data.arginfo import (
    init_bulk_mdata_arginfo,
    init_surf_mdata_arginfo,
    init_reaction_mdata_arginfo,
)

simplify_mdata_doc = simplify_mdata_arginfo().gen_doc(make_anchor=True)
with open('simplify-mdata-auto.rst', 'w') as f:
    f.write(simplify_mdata_doc)

init_bulk_mdata_doc = init_bulk_mdata_arginfo().gen_doc(make_anchor=True)
with open('init_bulk-mdata-auto.rst', 'w') as f:
    f.write(init_bulk_mdata_doc)

init_surf_mdata_doc = init_surf_mdata_arginfo().gen_doc(make_anchor=True)
with open('init_surf-mdata-auto.rst', 'w') as f:
    f.write(init_surf_mdata_doc)

init_reaction_mdata_doc = init_reaction_mdata_arginfo().gen_doc(make_anchor=True)
with open('init_reaction-mdata-auto.rst', 'w') as f:
    f.write(init_reaction_mdata_doc)

run_mdata_doc = run_mdata_arginfo().gen_doc(make_anchor=True)
with open('run-mdata-auto.rst', 'w') as f:
    f.write(run_mdata_doc)
