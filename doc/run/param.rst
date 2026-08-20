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

DPA4 and the PyTorch-exportable backend require DeePMD-kit 3.2 or later. A
regular PyTorch DPA4 model can be exported as an AOTInductor archive with:

.. code-block:: json

   {
     "train_backend": "pytorch",
     "model_format": "pt2"
   }

For DPA4C, and for DPA4 using the PyTorch-exportable backend, select the graph
export explicitly:

.. code-block:: json

   {
     "train_backend": "pytorch-exportable",
     "model_format": "pt2",
     "dp_compress": true
   }

``pt-expt`` is accepted as an alias of ``pytorch-exportable``. The ``pte``
format remains available for dense PyTorch-exportable models. Training
checkpoints keep the ``.pt`` suffix independently of the frozen model format.
