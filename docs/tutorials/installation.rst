Installation
============
Here are some advanced installation techniques, such as installing into a virtual environment.

.. warning::

	This program was developed for (and hopefully works on) \*nix and Mac computers.
	I've tried to make it work with Windows but haven't tested it; if you're using Windows
	and run into problems, please
	`open an issue on GitHub <https://github.com/adamjorr/kbbq/issues/new>`_.

Virtual Environments
--------------------

Virtual environments are separate environments that help prevent environment variables,
programs, and libraries from interfering with one another. An environment will have access
to system-wide programs but will use its own Python environment.
Installing into a virtual environment can help prevent a lot of future headaches and issues
that can be difficult to debug.

This tutorial will cover several methods for creating a virtual environment and installing
kbbq in the new virtual environment. I recommend using :ref:`inst_pipenv`, but the choice is yours.

 - :ref:`inst_pipenv`
 - :ref:`conda`
 - :ref:`venv`

.. _inst_pipenv:

Creating a Virtual Environment with pipenv
------------------------------------------

`Pipenv <https://docs.pipenv.org/>`_ is a venv-compatible virtual environment manager.
It's installable with pip::

	pip install --user pipenv

This will do a user install of pipenv, installing it in your home directory to prevent it
from interfering with system packages. Since it's a user install, you'll want to also add
the ``~/.local/bin`` directory to your ``PATH``. If you use ``bash`` as your shell, this can
be done by editing your ``~/.bash_profile``, then running::

	. ~/.bash_profile

To test that your installation worked properly and pipenv was added to your path, try
running ``pipenv --help``.

Now, change into the directory you want the virtual environment to live in; for example::

	mkdir myproject
	cd myproject

and install kbbq from git::

	pipenv install git+https://github.com/adamjorr/kbbq.git#egg=kbbq

.. note::

	We have to use the ``#egg=kbbq`` at the end of the git link because kbbq doesn't yet have
	a stable version number and hasn't been uploaded to pypi. Once it is, it will be
	installable with the command ``pipenv install kbbq``.

To work in this new environment, activate the shell with::

	pipenv shell

Now you're set! Type ``exit`` when you're done working to leave the virtual environment.
When you need to get back into the virtual environment, run the ``pipenv shell`` command
once again.

.. _conda:

Creating a Virtual Environment with conda
-----------------------------------------

.. note::

	I currently haven't written a conda package for kbbq. Once I do, instructions to install it
	will live here :)

.. _venv:

Creating a Virtual Environment with venv
----------------------------------------

:mod:`venv` is a python standard library module for creating virtual environments.
It's a low-level but lightweight virtual environment manager that comes with Python.

To make a virtual environment, enter the directory you want the virtual
environment to live in; for example::

	mkdir myproject
	cd myproject

Now pick a name for the environment and create it; here we will use the name ``env``::

	python3 -m venv env

This command tells python to use the main script of the :mod:`venv` module, and tells
:mod:`venv` that we want to make a virtual environment named env in the current directory.

To use the environment, you activate it by running the ``source`` command on the
``activate`` script :mod:`venv` created::

	source env/bin/activate

Once you do this, you should see the name of your environment on your command prompt.
You can now install any packages you want with pip and they will be installed into the
environment. To install kbbq::

	pip install git+https://github.com/adamjorr/kbbq.git

When you're finished using the environment, deactivate it by running::

	deactivate

Remember to reactivate your environment before doing more work on your project.

Testing that installation was successful
-----------------------------------------

To test that your installation was successful, run ::

	kbbq --help

while in your activated environment.
If you get usage information, installation probably succeeded.
