from dpgen.generator.arginfo import (
    run_mdata_arginfo,
    run_jdata_arginfo,
)

run_mdata_doc = run_mdata_arginfo().gen_doc(make_anchor=True)
with open('run-mdata-auto.rst', 'w') as f:
    f.write(run_mdata_doc)

run_jdata_doc = run_jdata_arginfo().gen_doc(make_anchor=True)
with open('run-jdata-auto.rst', 'w') as f:
    f.write(run_jdata_doc)
