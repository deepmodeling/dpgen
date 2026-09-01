=============================
dpgen run param parameters
=============================

.. note::
   One can load, modify, and export the input file by using our effective web-based tool `DP-GUI <https://dpgui.deepmodeling.com/input/dpgen-run>`_ online or hosted using the :ref:`command line interface <cli>` :code:`dpgen gui`. All parameters below can be set in DP-GUI. By clicking "SAVE JSON", one can download the input file.

.. dargs::
   :module: dpgen.generator.arginfo
   :func: run_jdata_arginfo

DPA4 and DPA4C
---------------

DPA4 and the PyTorch-exportable backend require DeePMD-kit 3.2 or later. DPA4
uses the regular PyTorch backend for both training and export:

.. code-block:: json

   {
     "train_backend": "pytorch",
     "model_format": "pt2",
     "default_training_param": {
       "model": {
         "type": "dpa4",
         "use_compile": true,
         "enable_tf32": true
       }
     }
   }

DPA4C uses the PyTorch-exportable backend for both training and graph export:

.. code-block:: json

   {
     "train_backend": "pytorch-exportable",
     "model_format": "pt2",
     "dp_compress": true,
     "default_training_param": {
       "model": {
         "descriptor": {"type": "dpa4c"}
       },
       "training": {
         "enable_compile": true,
         "enable_tf32": true
       }
     }
   }

The default ``train_backend`` remains ``tensorflow``. ``pt-expt`` is accepted
as an alias of ``pytorch-exportable``. PyTorch-exportable model deviation with
LAMMPS defaults to ``pt2``. Training checkpoints keep the ``.pt`` suffix
independently of the frozen model format.

Freeze and export use the same backend as training. Regular PyTorch and
PyTorch-exportable checkpoints are backend-specific and cannot be converted by
switching the ``dp`` backend flag after training.

The acceleration controls belong to different sections of the DeePMD training
template: DPA4 uses ``model.use_compile`` and ``model.enable_tf32``; DPA4C uses
``training.enable_compile`` and ``training.enable_tf32``. DP-GEN validates the
backend and these locations but does not inject either policy. A template cannot
mix DPA4 and DPA4C branches because they require different training backends.

AOTInductor ``.pt2`` archives are specific to the target CPU or GPU, GPU
compute capability, and libtorch version. DP-GEN therefore finishes the
training submission with the checkpoint, then runs ``freeze`` and optional
``compress`` in a separate submission using ``model_devi_machine`` and
``model_devi_resources``. Those resources must select the same hardware and
software target used by all subsequent model-deviation jobs. The resulting
models are linked into the model-deviation stage automatically.

The dense PyTorch-exportable ``pte`` format remains available for non-LAMMPS
workflows but is not supported by LAMMPS model deviation.
