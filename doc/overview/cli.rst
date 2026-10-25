.. _cli:

Command line interface
======================

.. argparse::
   :module: dpgen.main
   :func: main_parser
   :prog: dpgen

Offline workflow guides
-----------------------

Public workflow instructions from the repository's ``skills/`` directory are
bundled with installed DP-GEN releases. Read the guides for the installed
version without installing skills separately:

.. code-block:: console

   dpgen skills list
   dpgen skills read dpgen-run
   dpgen skills read dpgen-run references/monitoring.md
   dpgen skills read dpgen-simplify assets/param.template.json

``dpgen skills read`` accepts an optional relative text-file path within the
selected guide, defaulting to ``SKILL.md``. These commands work offline and
only display documentation: they do not execute examples, start workflows,
or submit HPC jobs.
