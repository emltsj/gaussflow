# Contributing to GaussFlow

Thanks for your interest in improving GaussFlow.

## Development setup

```bash
git clone https://github.com/emltsj/gaussflow.git
cd gaussflow
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python -m unittest discover -s tests -v
```

## Guidelines

- Keep the deterministic core deterministic: the LLM may propose structures or
  calculation chains, but input writing, parsing, unit conversion, validation,
  and export must remain program-controlled.
- Do not commit paper-specific production protocols, local absolute paths,
  API keys, Gaussian outputs, run directories, databases, or logs.
- Add or update tests for parser, validation, recipe, or recovery changes.
- Keep new calculation chains generic and configurable (methods, basis sets,
  solvents, stages).
- Run the full test suite before opening a pull request.

## Reporting bugs

Please include:

- GaussFlow version / commit
- Windows version and Python version
- The command or GUI action
- The relevant `gflow ... --json` output, without API keys or private structures
- A short excerpt of the Gaussian output if relevant

Do not paste API keys, personal paths, or unpublished molecular data into issues.

## License

By contributing, you agree that your contributions are licensed under the MIT
License (see [LICENSE](LICENSE)).
