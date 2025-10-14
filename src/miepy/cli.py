"""CLI entry point for running miepy examples."""

import sys
from pathlib import Path


def get_examples_dir() -> Path:
    """Get the examples directory, works both in dev and installed package."""
    # Examples are now part of the package
    examples_dir = Path(__file__).parent / "examples"
    if examples_dir.exists():
        return examples_dir

    raise FileNotFoundError("Could not find examples directory")


def main() -> None:
    """Run a miepy example script by name."""
    if len(sys.argv) < 2:
        print("Usage: uvx miepy <example-name>")
        print("\nAvailable examples:")
        try:
            examples_dir = get_examples_dir()
            examples = {}
            # Map user-friendly names to file names
            for ex in sorted(examples_dir.glob("*.py")):
                if ex.stem != "__init__":
                    # Remove number prefix for display
                    name = ex.stem.lstrip("0123456789_")
                    examples[name] = ex

            for name in sorted(examples.keys()):
                print(f"  - {name}")
        except FileNotFoundError:
            print("  (examples not found)")
        sys.exit(1)

    example_name = sys.argv[1]

    try:
        examples_dir = get_examples_dir()
    except FileNotFoundError:
        print("Error: Examples directory not found")
        sys.exit(1)

    # Try to find example by name (with or without number prefix)
    example_file = None
    for ex in examples_dir.glob("*.py"):
        if ex.stem != "__init__":
            name = ex.stem.lstrip("0123456789_")
            if name == example_name or ex.stem == example_name:
                example_file = ex
                break

    if example_file is None or not example_file.exists():
        print(f"Error: Example '{example_name}' not found")
        sys.exit(1)

    # Execute the example file directly
    with example_file.open() as f:
        code = compile(f.read(), str(example_file), "exec")
        exec(code, {"__name__": "__main__"})


if __name__ == "__main__":
    main()
