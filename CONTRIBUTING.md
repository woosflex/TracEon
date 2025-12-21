# Contributing to TracEon

We welcome contributions! Please follow these guidelines to ensure a smooth process.

## How to Contribute
1. **Fork** the repository.
2. **Create a branch** for your feature or bug fix.
3. **Commit** your changes with clear messages.
4. **Push** to your fork and submit a **Pull Request (PR)**.

## Requirements
- **Tests:** New features must include unit tests. Run `ctest` to verify.
- **Style:** Follow modern C++20 conventions (RAII, `std::string_view`, etc.).
- **Performance:** If you modify core logic, run the regression check:
  ```bash
  python benchmarks/check_regression.py --baseline main
  ```

## Reporting Issues
Please include:
- **System details** (OS, CPU, Compiler).
- A minimal **reproducible example**.