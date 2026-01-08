# Contributing to CAL-DISP

Thanks for contributing! This guide will help you get started.

## Development Setup
```bash
# Clone and install
git clone https://github.com/opera-adt/cal-disp.git
cd cal-disp
mamba env create -f environment.yml
conda activate my-cal-env

# Install with dev dependencies
python -m pip install -e ".[download,test]"

# Set up pre-commit hooks
pre-commit install
```

## Code Standards

We keep things simple and consistent:

**Style:**
- Use [black](https://black.readthedocs.io/) for formatting (88 char line length)
- Use [ruff](https://docs.astral.sh/ruff/) for linting
- Use [mypy](https://mypy-lang.org/) for type checking
- Sort imports with [isort](https://pycqa.github.io/isort/)

**Docstrings:**
- Follow [numpydoc format](https://numpydoc.readthedocs.io/en/latest/format.html)
- Include Parameters, Returns, Raises sections
- Add examples for public functions

**Type Hints:**
- Required for all public functions
- Use `from __future__ import annotations` for forward references
- Example: `def process(data: Path, output: Path | None = None) -> None:`

**Code Style:**
- Write clear, readable code over clever code
- Keep functions focused and short
- Prefer explicit over implicit
- Use descriptive names for broad scopes, short names for tight loops

## Making Changes

1. **Create a branch:**
```bash
   git checkout -b feature/my-feature
```

2. **Make your changes:**
   - Write tests for new functionality
   - Update docstrings
   - Keep commits focused and atomic

3. **Test your changes:**
```bash
   # Run tests
   pytest

   # Run specific test file
   pytest tests/test_myfile.py

   # Run with coverage
   pytest --cov=cal_disp --cov-report=html

   # Check types
   mypy src/
```

4. **Commit:**
```bash
   git add .
   git commit -m "Add feature: brief description"
```

   Pre-commit hooks will automatically run. If they fail:
   - Fix the issues
   - Re-add modified files: `git add .`
   - Commit again

5. **Push and create PR:**
```bash
   git push origin feature/my-feature
```
   Then open a pull request on GitHub.

## Pull Request Guidelines

**Before submitting:**
- [ ] Tests pass locally (`pytest`)
- [ ] Pre-commit checks pass
- [ ] Code is documented
- [ ] Changes are described in PR description

**PR description should include:**
- What problem does this solve?
- How does it solve it?
- Any breaking changes?
- Related issues (use `Fixes #123` to auto-close)

## Running Tests
```bash
# All tests
pytest

# Specific test
pytest tests/test_download.py::test_frame_download

# With verbose output
pytest -v

# Stop on first failure
pytest --maxfail=1

# Run only tests matching pattern
pytest -k "download"
```

## Common Tasks

**Update dependencies:**
```bash
# Update environment.yml, then:
mamba env update -f environment.yml --prune
```

**Add new dependency:**
1. Add to `environment.yml` (conda) or `pyproject.toml` (pip)
2. Reinstall: `mamba env update -f environment.yml --prune`

**Debug test failures:**
```bash
# Run with print statements visible
pytest -s

# Drop into debugger on failure
pytest --pdb
```

## Code Review Process

1. Maintainers review PRs within a few days
2. Address feedback by pushing new commits
3. Once approved, maintainers will merge
4. PRs require passing CI and at least one approval

## Questions?

- Open an issue for bugs or feature requests
- Tag maintainers in PR comments for urgent questions
- Check existing issues and PRs first

## License

By contributing, you agree that your contributions will be licensed under the same license as the project (BSD-3-Clause or Apache-2.0).
