Development Guide
==================

Found a bug? Saw a typo? Want a new feature? Any questions or suggestions? Please don't be shy! Tell us by opening an issue or just start a pull request.

If you want to contribute to spatialtis, here is some information to start with.

First clone the repository::

    git clone https://github.com/Mr-Milk/SpatialTis.git

And then install the dependencies, it should contain everything::

    pip install -r requirements.txt

Here are some tools used in development:

    Code Quality:

    - `pre-commit hook <https://pre-commit.com/>`_
        - `isort <https://pycqa.github.io/isort/>`_ (sort import)
        - `black <https://black.readthedocs.io/en/stable/>`_ (code formating)
        - `flake8 <https://flake8.pycqa.org/en/latest/>`_ (pylint check)
    - `pytest <https://docs.pytest.org/en/latest/>`_ (test framework)

    Documentation:

    - `sphinx <https://www.sphinx-doc.org/en/master/>`_

    Continuous Integration:

    - `Github Actions <https://github.com/features/actions>`_

And here are some requirements:

    - Ideally, the test coverage should be > 90%.
    - Every exposed API should have type annotation.
    - Docstring is google style

