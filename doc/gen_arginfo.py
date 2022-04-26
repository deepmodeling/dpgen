from dpgen.generator.arginfo import run_mdata_arginfo

run_mdata_doc = run_mdata_arginfo().gen_doc()
with open('run-mdata-auto.rst', 'w') as f:
    f.write(run_mdata_doc)
