# Contributing

## Development setup

```bash
git clone https://github.com/imuhdawood/scinterop.git
cd scinterop
pip install -e .[dev]
```

## Running tests

```bash
python -m pytest
```

## Building docs

```bash
pip install mkdocs mkdocs-material mkdocstrings[python]
mkdocs serve
```

## Code style

- Google-style docstrings
- Follow existing patterns

## Submitting a PR

1. Fork the repo.
2. Create a feature branch.
3. Add tests for any new functionality.
4. Run all tests with `python -m pytest`.
5. Submit a PR against `main`.
