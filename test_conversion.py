"""
Test script to verify the Python conversion structure.

This script validates the code structure without requiring data files
or heavy dependencies to be installed.
"""

import ast
import sys


def check_python_file(filepath):
    """Check if a Python file is syntactically valid."""
    try:
        with open(filepath, 'r') as f:
            code = f.read()
        ast.parse(code)
        print(f"✓ {filepath}: Syntax valid")
        return True
    except SyntaxError as e:
        print(f"✗ {filepath}: Syntax error at line {e.lineno}: {e.msg}")
        return False


def check_imports(filepath):
    """Extract and display imports from a Python file."""
    with open(filepath, 'r') as f:
        code = f.read()
    
    tree = ast.parse(code)
    imports = []
    
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imports.append(alias.name)
        elif isinstance(node, ast.ImportFrom):
            # ImportFrom can have None module for relative imports
            if node.module:
                imports.append(node.module)
    
    return list(set(imports))


def check_functions(filepath):
    """Extract function names from a Python file."""
    with open(filepath, 'r') as f:
        code = f.read()
    
    tree = ast.parse(code)
    functions = []
    
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            functions.append(node.name)
    
    return functions


def main():
    print("=" * 70)
    print("Python Conversion Validation")
    print("=" * 70)
    
    files_to_check = [
        'geological_analysis.py',
        'visualization_utils.py'
    ]
    
    all_valid = True
    
    for filepath in files_to_check:
        print(f"\n{'=' * 70}")
        print(f"Checking: {filepath}")
        print('=' * 70)
        
        # Check syntax
        if not check_python_file(filepath):
            all_valid = False
            continue
        
        # Check imports
        imports = check_imports(filepath)
        print(f"\nImported modules ({len(imports)}):")
        for imp in sorted(imports):
            if imp:
                print(f"  - {imp}")
        
        # Check functions
        functions = check_functions(filepath)
        print(f"\nDefined functions ({len(functions)}):")
        for func in functions[:10]:  # Show first 10
            print(f"  - {func}()")
        if len(functions) > 10:
            print(f"  ... and {len(functions) - 10} more")
    
    print(f"\n{'=' * 70}")
    if all_valid:
        print("✓ All files are syntactically valid!")
        print("✓ Python conversion structure is correct")
        print("\nNext steps:")
        print("1. Install dependencies: pip install -r requirements.txt")
        print("2. Add data files to data/ directory")
        print("3. Run: python geological_analysis.py")
    else:
        print("✗ Some files have syntax errors")
        sys.exit(1)
    print("=" * 70)


if __name__ == "__main__":
    main()
