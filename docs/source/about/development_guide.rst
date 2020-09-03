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
    - pre-commit framework

        - isort (sort import)
        - black (code formating)
        - mypy (type check)
        - flake8 (pylint check)

    - pytest (test framework)

    Documentation:
    - sphinx

    Continuous Integration:
    - Github Actions

And here are some requirements:

    - Ideally, the test coverage should be > 80%.
    - Every exposed API should have type annotation.
    - Docstring is google style

